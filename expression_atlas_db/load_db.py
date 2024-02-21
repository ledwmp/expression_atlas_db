import sys
from pathlib import Path
from typing import Union
import click
import logging
import time

def configure_logger(log_file=None, level=logging.INFO, overwrite_log=True,
                     format=logging.BASIC_FORMAT):
    # console and file
    if log_file is None:
        logging.basicConfig(stream=sys.stdout, level=level, format=format)
    else:
        logging.basicConfig(filename=log_file, level=level,
                            filemode=('w' if overwrite_log else 'a'),
                            format=format)
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(format))
        logging.getLogger('').addHandler(console)
configure_logger(f"{time.strftime('%Y%m%d_%H%M%S')}_veliadb_load_db.log",
                 level=logging.INFO)


from expression_atlas_db import base, settings
from expression_atlas_db.utils import GTFParser, ExperimentParser, MetaDataFetcher

def bulk_insert_gtf(session:base._Session, gtf:Union[GTFParser,None], batch_columns:int=10000) -> None:
    """Read from gtf and populate the transcript and gene tables in expression_atlas_db.
    gtf should be an object with transcript_df and gene_df.
    TODO: Replace None above with other sources, perhaps velia_db.

    Args:
        session (base._Session) base session for expression_atlas_db.
        gtf (Union[GTFParser,None]) object holding transcript_df and gene_df.
        batch_columns (int) number of columns to batch at one time when running inserts.
    """
    # Add all genes from veliadb into sequenceregion, then gene.

    logging.info('Creating gene sequenceregions.')
    records = [f'({r.veliadb_id}, "{r.assembly_id}", "{"gene"}")' for i,r in gtf.gene_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'sequenceregion' ('veliadb_id','assembly_id','type') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

    gene_srs = { g.veliadb_id:g for g in session \
                    .query(base.SequenceRegion) \
                    .filter((base.SequenceRegion.type == 'gene')).all()}
    
    logging.info('Populating gene table.')
    records = [f'({gene_srs[r.veliadb_id].id}, "{r.gene_id}", "{r.gene_biotype}")' for i,r in gtf.gene_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'gene' ('id','gene_id','gene_biotype') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};""")
        session.commit()

    genes = {g.gene_id:g for g in session.query(base.Gene).all()}
    
    # Add all transcripts from veliadb into sequenceregion, then transcript.

    logging.info('Creating transcript sequence regions.')
    records = [f'({r.veliadb_id}, "{r.assembly_id}", "{"transcript"}")' for i,r in gtf.transcript_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'sequenceregion' ('veliadb_id','assembly_id','type') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

    transcript_srs = {t.veliadb_id:t for t in session \
                        .query(base.SequenceRegion) \
                        .filter((base.SequenceRegion.type == 'transcript')).all()}

    logging.info('Populating transcript table.')
    records = [f'({transcript_srs[r.veliadb_id].id}, "{r.transcript_id}", {genes[r.gene_id].id}, "{r.gene_biotype}")' for i,r in gtf.transcript_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'transcript' ('id','transcript_id','gene_id','gene_biotype') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
            )
        session.commit()

def bulk_insert_measurements(session:base._Session, meta:MetaDataFetcher, exp:ExperimentParser, batch_columns:int=10000,
) -> None:
    """Bulk insert studies into expression_atlas_db. Some of the functionality here assumes single-threaded 
    access to the database and lack of any write/deletes during the load. 

    Args:  
        session (base._Session)
        meta (MetaDataFetcher)
        exp (ExperimentParser)
    """
    
    logging.info(f'Adding study {meta._study_id}.')
    study = base.Study(
                geo_id=meta._geo_id,
                srp_id=meta._srp_id, 
                pmid=','.join(meta._pmids), 
                bio_id=meta._project_id,
                title=meta._project_title,
                description=meta._project_summary,
                )
    
    session.add(study)
    session.commit()

    samples_metadata = meta.samples_metadata
    samples_metadata = samples_metadata.merge(
        exp.samples_metadata[exp.samples_metadata.columns[~exp.samples_metadata.columns.isin(samples_metadata.columns)]],
        left_on='Experiment',
        right_index=True,
        )
    
    logging.info(f'Adding {samples_metadata.shape[0]} samples from {meta._study_id}.')
    samples = [base.Sample(
                        srx_id=s['Experiment'], 
                        atlas_group=s['sample_type_1'], 
                        study=study, 
                        fields=s,
                        ) for s in samples_metadata.to_dict('records')]
    session.add_all(samples)
    session.commit()

    sequenceregions = {g.gene_id:g for g in session.query(base.Gene).all()}
    transcripts = {t.transcript_id:t for t in session.query(base.Transcript).all()}
    sequenceregions.update(transcripts)

    gm_samples, gm_regions, gm_measurments, measurements = exp.prepare_measurements(measurement_type='gene')
    g_exp_dict = exp.prepare_differentialexpression(measurement_type='gene')
    tm_samples, tm_regions, tm_measurments, measurements = exp.prepare_measurements(measurement_type='transcript')
    t_exp_dict = exp.prepare_differentialexpression(measurement_type='transcript')

    measurement_convert = {
                    'counts':'counts',
                    'normed_counts':'normed_counts',
                    'raw_tpm':'tpm',
                    'normed_counts_transform':'normed_counts_transform'
                    }
    measurement_columns = [measurement_convert[m] for m in measurements]
    de_columns = ['basemean', 'log2foldchange', 'lfcse', 'stat', 'pvalue', 'padj', 'log10_padj', 'control_mean', 'case_mean']

    for c, (c_df, c_left, c_right) in g_exp_dict.items():

        logging.info('Adding contrast {c} from study {study.srp_id}.')

        left_condition_display = c_df['left_condition_display'].unique()[0]
        right_condition_display = c_df['right_condition_display'].unique()[0]
        sample_condition_key = c_df['sample_condition_key'].unique()[0]
        left_condition_level = c_df['left_condition_level'].unique()[0]
        right_condition_level = c_df['right_condition_level'].unique()[0]

        contrast = base.Contrast(
                        study=study, 
                        left_condition_display=left_condition_display, 
                        right_condition_display=right_condition_display,
                        sample_condition_key=sample_condition_key,
                        left_condition_level=left_condition_level,
                        right_condition_level=right_condition_level,
                        contrast_name=c,
                        )
        session.add(contrast)
        session.commit()

        logging.info('Adding left samplecontrasts from study {study.srp_id}.')
        left_samples = session.query(base.Sample).filter(base.Sample.srx_id.in_(c_left)).all()
        left_sample_contrasts = [base.SampleContrast(
                                                sample=s, 
                                                contrast_side='left', 
                                                contrast=contrast, 
                                                study=study,
                                                condition_display=contrast.left_condition_display,
                                                condition_level=contrast.left_condition_level,
                                                sample_condition_key=contrast.sample_condition_key,
                                                ) for s in left_samples]

        logging.info('Adding right samplecontrasts from study {study.srp_id}.')
        right_samples = session.query(base.Sample).filter(base.Sample.srx_id.in_(c_right)).all()
        right_sample_contrasts = [base.SampleContrast(
                                                sample=s, 
                                                contrast_side='right', 
                                                contrast=contrast, 
                                                study=study,
                                                condition_display=contrast.right_condition_display,
                                                condition_level=contrast.right_condition_level,
                                                sample_condition_key=contrast.sample_condition_key,
                                                ) for s in right_samples]
        session.add_all(left_sample_contrasts+right_sample_contrasts)
        session.commit()

        logging.info('Adding measurements for gene de from study {study.srp_id}.')
        records = [f'({sequenceregions[i].id}, "{study.id}", "{"differentialexpression"}")' for i,r in c_df.iterrows()]
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                    f"""INSERT INTO 'measurement' ('sequenceregion_id','study_id','type') VALUES 
                    {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                    )
            session.commit()

        de_ms = [de.id for de in session.query(base.Measurement).filter( \
                        (base.Measurement.type == 'differentialexpression') & \
                        (base.Measurement.study_id == study.id)
                        ).order_by(base.Measurement.id.desc()).limit(c_df.shape[0])
                    ][::-1]

        records = [f'({k}, {contrast.id}, {", ".join(map(lambda x: str(x),[r[c] for c in de_columns]))})' \
                    for (i,r),k in zip(c_df.iterrows(), de_ms)]

        logging.info('Adding gene de measurements from study {study.srp_id}.')
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                    f"""INSERT INTO 'differentialexpression' ("id", "contrast_id", {', '.join(map(lambda x: "'"+x+"'", de_columns))}) VALUES 
                    {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                    )
            session.commit()

    logging.info('Adding measurements for gene samples from study {study.srp_id}.')
    records = [f'({sequenceregions[i].id}, "{study.id}", "{"samplemeasurement"}")' for i in gm_regions]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'measurement' ('sequenceregion_id','study_id','type') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

    sample_ms = [sm.id for sm in session.query(base.Measurement).filter( \
                    (base.Measurement.type == 'samplemeasurement') & \
                    (base.Measurement.study_id == study.id)
                    ).order_by(base.Measurement.id.desc()).limit(gm_regions.shape[0])
                ][::-1]

    samples = {s.srx_id:s for s in session.query(base.Sample).filter(base.Sample.srx_id.in_(gm_samples))}

    records = [f'({i}, {samples[s].id}, {", ".join(map(lambda x: str(x), a))})' \
                            for i,s,a in zip(sample_ms, gm_samples, gm_measurments)]

    logging.info('Adding gene sample measurments from study {study.srp_id}.')
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'samplemeasurement' ("id", "sample_id", {', '.join(map(lambda x: "'"+x+"'", measurement_columns))}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

    for c, (c_df, c_left, c_right) in t_exp_dict.items():

        logging.info('Adding measurements for transcript de from study {study.srp_id}.')
        records = [f'({sequenceregions[i].id}, "{study.id}", "{"differentialexpression"}")' for i,r in c_df.iterrows()]
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                    f"""INSERT INTO 'measurement' ('sequenceregion_id','study_id','type') VALUES 
                    {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                    )
            session.commit()

        de_ms = [de.id for de in session.query(base.Measurement).filter( \
                        (base.Measurement.type == 'differentialexpression') & \
                        (base.Measurement.study_id == study.id)
                        ).order_by(base.Measurement.id.desc()).limit(c_df.shape[0])
                    ][::-1]

        records = [f'({k}, {contrast.id}, {", ".join(map(lambda x: str(x),[r[c] for c in de_columns]))})' \
                    for (i,r),k in zip(c_df.iterrows(), de_ms)]

        logging.info('Adding transcript de measurements from study {study.srp_id}.')
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                    f"""INSERT INTO 'differentialexpression' ("id", "contrast_id", {', '.join(map(lambda x: "'"+x+"'", de_columns))}) VALUES 
                    {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                    )
            session.commit()

    logging.info('Adding measurements for transcript samples from study {study.srp_id}.')
    records = [f'({sequenceregions[i].id}, "{study.id}", "{"samplemeasurement"}")' for i in tm_regions]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'measurement' ('sequenceregion_id','study_id','type') VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

    sample_ms = [sm.id for sm in session.query(base.Measurement).filter( \
                    (base.Measurement.type == 'samplemeasurement') & \
                    (base.Measurement.study_id == study.id)
                    ).order_by(base.Measurement.id.desc()).limit(tm_regions.shape[0])
                ][::-1]

    samples = {s.srx_id:s for s in session.query(base.Sample).filter(base.Sample.srx_id.in_(tm_samples))}

    records = [f'({i}, {samples[s].id}, {", ".join(map(lambda x: str(x), a))})' \
                            for i,s,a in zip(sample_ms, tm_samples, tm_measurments)]

    logging.info('Adding transcript sample measurements from study {study.srp_id}.')
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO 'samplemeasurement' ("id", "sample_id", {', '.join(map(lambda x: "'"+x+"'", measurement_columns))}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
        session.commit()

@click.command()
@click.option('--drop-all', is_flag=True, help='Empty database and reload data.')
def load_db(drop_all):
    """
    """
    if drop_all:
        logging.info('Dropping the database')
        base.Base.metadata.drop_all()
        
        logging.info('Creating the database')
        base.Base.metadata.create_all()

    # make the session
    session = base.Session()


if __name__ == '__main__':
    load_db()