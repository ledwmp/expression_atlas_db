import sys
from pathlib import Path
from typing import Union
import click
import logging
import time
import numpy as np
import pandas as pd
from typing import Dict, List, Union
import pathlib
import s3fs

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

from expression_atlas_db import base, settings, redshift_stmts
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
    records = [f"({r.veliadb_id}, '{r.assembly_id}', '{'gene'}')" for i,r in gtf.gene_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO sequenceregion (veliadb_id, assembly_id, type) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )

    gene_srs = { g.veliadb_id:g for g in session \
                    .query(base.SequenceRegion) \
                    .filter((base.SequenceRegion.type == 'gene')).all()}
    
    logging.info('Populating gene table.')
    records = [f"({gene_srs[r.veliadb_id].id}, '{r.gene_id}', '{r.gene_biotype}')" for i,r in gtf.gene_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO gene (id, gene_id, gene_biotype) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};""")

    genes = {g.gene_id:g for g in session.query(base.Gene).all()}
    
    # Add all transcripts from veliadb into sequenceregion, then transcript.

    logging.info('Creating transcript sequence regions.')
    records = [f"({r.veliadb_id}, '{r.assembly_id}', '{'transcript'}')" for i,r in gtf.transcript_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO sequenceregion (veliadb_id, assembly_id , type) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )

    transcript_srs = {t.veliadb_id:t for t in session \
                        .query(base.SequenceRegion) \
                        .filter((base.SequenceRegion.type == 'transcript')).all()}

    logging.info('Populating transcript table.')
    records = [f"({transcript_srs[r.veliadb_id].id}, '{r.transcript_id}', {genes[r.gene_id].id}, '{r.gene_biotype}')" for i,r in gtf.transcript_df.iterrows()]
    for i in range(0, int(len(records)/batch_columns)+1):
        session.execute(
                f"""INSERT INTO transcript (id, transcript_id, gene_id, gene_biotype) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
            )
    session.commit()

def create_differentialexpression(
                    de_df: pd.DataFrame,
                    sequenceregions: Dict[str,base.SequenceRegion],
                    study: base.Study,
                    contrast: base.Contrast,
                    batch_columns:int=10000,
                    bulk_insert:bool=False,
                    local_dump:bool=False,
                    s3_dump:bool=False,
                    session: Union[base._Session,None]=None,
                    de_columns:List[str]=[c.name for c in base.DifferentialExpression.__table__.columns \
                                                    if c.name not in ('id','contrast_id','sequenceregion_id')],
                    s3_fs:Union[s3fs.core.S3FileSystem,None]=None,
                    fh:Union[pathlib.PosixPath,None]=None,
                    ) -> None:
    """ 
    Args:
        de_df (pd.DataFrame)
        sequenceregions (Dict[str,base.SequenceRegion])
        study (base.Study)
        contrast (base.Contrast)
        batch_columns (int)
        bulk_insert (bool)
        local_dump (bool)
        s3_dump (bool)
        session (Union[base._Session,None])
        de_coulumns (List[str])
        s3_fs (s3fs.core.S3FileSystem)
        fh (Union[pathlib.PosixPath,None])
    """
    if bulk_insert and not session:
        raise ValueError('Need to provide a session if using bulk insert.')
    if s3_dump and not (s3_fs and fh):
        raise ValueError('Need to pass in s3fs object and s3 location if using s3 dump.')
    if local_dump and not fh:
        raise ValueError('Need to pass in fh location if using local dump.')
    if bulk_insert:
        records = [f'({contrast.id}, {sequenceregions[i].id}, {", ".join(map(lambda x: str(x),[r[c] for c in de_columns]))})' \
                                                                        for i,r in de_df.iterrows()]

        logging.info(f'Adding de measurements from study {study.srp_id} via bulk insert.')
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                f"""INSERT INTO differentialexpression (contrast_id, sequenceregion_id, {', '.join(de_columns)}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
            session.commit()
    else:
        records = [','.join(map(lambda x: str(x), [contrast.id, sequenceregions[i].id, *[r[c] for c in de_columns]])) for i, r in de_df.iterrows()]
        logging.info(f'Adding de measurments from study {study.srp_id} to csv at {fh}.')
        if local_dump:
            with open(fh, 'w') as f_out:
                f_out.write('\n'.join(records))
        elif s3_dump:
            # Check that s3 location exists and is accessible.
            if not s3_fs.exists(str(fh.parent).replace('s3:/','s3://')):
                raise FileExistsError(f'Location for s3_dump {fh.parent} does not exists.')
            with s3_fs.open(str(fh).replace('s3:/','s3://'), 'w') as f_out:
                f_out.write('\n'.join(records))

def create_samplemeasurements(
                    m_regions: np.ndarray,
                    m_samples: np.ndarray,
                    m_measurements: np.ndarray,
                    sequenceregions: Dict[str,base.SequenceRegion],
                    study: base.Study,
                    samples: Dict[str,base.Sample],
                    batch_columns:int=10000,
                    bulk_insert:bool=False,
                    local_dump:bool=False,
                    s3_dump:bool=False,
                    session: Union[base._Session,None]=None,
                    measurement_columns:List[str]=[c.name for c in base.SampleMeasurement.__table__.columns \
                                                    if c.name not in ('id','sample_id','sequenceregion_id')],
                    s3_fs:Union[s3fs.core.S3FileSystem,None]=None,
                    fh:Union[pathlib.PosixPath,None]=None,
                    ) -> None:
    """ 
    Args:
        m_regions (np.ndarray)
        m_samples (np.ndarray)
        m_measurements (np.ndarray)
        sequenceregions (Dict[str,base.SequenceRegion])
        study (base.Study)
        samples (Dict[str,base.Sample])
        batch_columns (int)
    """
    if bulk_insert and not session:
        raise ValueError('Need to provide a session if using bulk insert.')
    if s3_dump and not (s3_fs and fh):
        raise ValueError('Need to pass in s3fs object and s3 location if using s3 dump.')
    if local_dump and not fh:
        raise ValueError('Need to pass in fh location if using local dump.')
    if bulk_insert:
        records = [f'({samples[s].id}, {sequenceregions[i].id}, {", ".join(map(lambda x: str(x), a))})' \
                                                    for i,s,a in zip(m_regions, m_samples, m_measurements)]

        logging.info(f'Adding sample measurments from study {study.srp_id} via bulk insert.')
        for i in range(0, int(len(records)/batch_columns)+1):
            session.execute(
                f"""INSERT INTO samplemeasurement (sample_id, sequenceregion_id, {', '.join(measurement_columns)}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
                )
            session.commit()
    else:
        records = [','.join(map(lambda x: str(x), [samples[s].id, sequenceregions[i].id, *a])) for i,s,a in zip(m_regions, m_samples, m_measurements)]
        logging.info(f'Adding sample measurments from study {study.srp_id} to csv at {fh}.')
        if local_dump:
            with open(fh, 'w') as f_out:
                f_out.write('\n'.join(records))
        elif s3_dump:
            # Check that s3 location exists and is accessible.
            if not s3_fs.exists(str(fh.parent).replace('s3:/','s3://')):
                raise FileExistsError(f'Location for s3_dump {fh.parent} does not exists.')
            with s3_fs.open(str(fh).replace('s3:/','s3://'), 'w') as f_out:
                f_out.write('\n'.join(records))

def copy_into_redshift_table(
                            session: base._Session,
                            table: base.Base, 
                            s3_staging_loc: pathlib.PosixPath,
                            iam_role:str=settings.redshift_iam_role,
                            ) -> None:
    """ 
    Args:
        session (base._Session) sqlalchemy session.
        table (base.Base) table class for destination table.
        s3_staging_loc (pathlib.PosixPath) the location of the staging directory where csv tables get copied.
        iam_role (str) iam_role string for redshift.
    """
    s3_copy_command = f"""copy {table.__tablename__} ({', '.join([c.name for c in table.__table__.columns if not c.primary_key])})
    from '{str(s3_staging_loc).replace('s3:/','s3://')}'
    iam_role '{iam_role}'
    delimiter ',';
    """
    try:
        logging.info(f'Copying table {s3_staging_loc} into {table.__tablename__}.')
        session.execute(s3_copy_command)
        session.commit()
    except Exception as e:
        logging.exception(e)

def bulk_insert_measurements(
                        session: base._Session, 
                        meta: MetaDataFetcher, 
                        exp: ExperimentParser,
                        batch_columns: int=100000,
                        s3: Union[s3fs.core.S3FileSystem,None]=None,
                        staging_loc: pathlib.PosixPath=settings.test_staging_loc,
                        ) -> None:
    """ 
    Args:
        session (base._Session)
        meta (MetaDataFetcher)
        exp (ExperimentParser)
        batch_columns (int)
        
    """

    logging.info(f'Adding study {meta._study_id}.')
    study = base.Study(
                velia_id=meta._study_id,
                geo_id=meta._geo_id,
                srp_id=meta._srp_id, 
                pmid=','.join(meta._pmids), 
                bio_id=meta._project_id,
                timestamps=exp.file_timestamps,
                sizes=exp.file_sizes,
                title=meta._project_title,
                description=meta._project_summary,
                )
    
    session.add(study)

    samples_metadata = meta.samples_metadata.merge(
        exp.samples_metadata[exp.samples_metadata.columns[~exp.samples_metadata.columns.isin(meta.samples_metadata.columns)]],
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

    sequenceregions = {g.gene_id:g for g in session.query(base.Gene).all()}
    transcripts = {t.transcript_id:t for t in session.query(base.Transcript).all()}
    sequenceregions.update(transcripts)

    measurements = [c.name for c in base.SampleMeasurement.__table__.columns \
                                                    if c.name not in ('id','sample_id','sequenceregion_id')]

    measurement_convert = {
                    'counts':'counts',
                    'normed_counts':'normed_counts',
                    'tpm':'raw_tpm',
                    'normed_counts_transform':'normed_counts_transform'
                    }
    measurement_columns = [measurement_convert[m] for m in measurements]

    gm_samples, gm_regions, gm_measurments, g_measurements = exp.prepare_measurements(measurement_type='gene', measurements=measurement_columns)
    g_exp_dict = exp.prepare_differentialexpression(measurement_type='gene')
    tm_samples, tm_regions, tm_measurments, t_measurements = exp.prepare_measurements(measurement_type='transcript', measurements=measurement_columns)
    t_exp_dict = exp.prepare_differentialexpression(measurement_type='transcript')

    for c, (c_df, c_left, c_right) in g_exp_dict.items():

        logging.info(f'Adding contrast {c} from study {study.srp_id}.')

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

        logging.info(f'Adding left samplecontrasts from study {study.srp_id}.')
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

        logging.info(f'Adding right samplecontrasts from study {study.srp_id}.')
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
        
        # create_differentialexpression(c_df, sequenceregions, study, contrast, local_dump=True, fh=f'{study.id}.{contrast.id}.gene_de.csv')
        create_differentialexpression(
                                    c_df, 
                                    sequenceregions, 
                                    study, 
                                    contrast, 
                                    s3_dump=True, 
                                    s3_fs=s3, 
                                    fh=staging_loc / f'gene_de.{study.id}.{contrast.id}.csv',
                                    )

    samples = {s.srx_id:s for s in session.query(base.Sample).filter(base.Sample.srx_id.in_(np.unique(gm_samples)))}
    # create_samplemeasurements(gm_regions, gm_samples, gm_measurments, sequenceregions, study, samples, local_dump=True, fh=f'{study.id}.gene_measurment.csv')
    create_samplemeasurements(
                            gm_regions, 
                            gm_samples, 
                            gm_measurments, 
                            sequenceregions, 
                            study, 
                            samples, 
                            s3_dump=True, 
                            s3_fs=s3,
                            fh=staging_loc / f'gene_measurement.{study.id}.csv',
                            )

    for c, (c_df, c_left, c_right) in t_exp_dict.items():
            
            contrast = session.query(base.Contrast).filter((base.Contrast.contrast_name == c) & (base.Contrast.study_id == study.id)).first()

            create_differentialexpression(
                                    c_df, 
                                    sequenceregions, 
                                    study, 
                                    contrast, 
                                    s3_dump=True, 
                                    local_dump=False,
                                    s3_fs=s3, 
                                    fh=staging_loc / f'transcript_de.{study.id}.{contrast.id}.csv',
                                    )
        
    samples = {s.srx_id:s for s in session.query(base.Sample).filter(base.Sample.srx_id.in_(np.unique(tm_samples)))}
    # create_samplemeasurements(tm_regions, tm_samples, tm_measurments, sequenceregions, study, samples, local_dump=True, fh=f'{study.id}.transcript_measurment.csv')
    create_samplemeasurements(
                            tm_regions, 
                            tm_samples, 
                            tm_measurments, 
                            sequenceregions, 
                            study, 
                            samples, 
                            s3_dump=True, 
                            local_dump=False,
                            s3_fs=s3,
                            fh=staging_loc / f'transcript_measurement.{study.id}.csv',
                            )
    session.commit()

@click.command()
@click.option('--drop-all', is_flag=True, help='Empty database and reload data.')
@click.option('--drop-dataset', is_flag=True, help='Drop dataset tablees.')
@click.option('-gtf', help='Path to gtf file.')
@click.option('--redshift', is_flag=True)
def load_db(gtf:str, drop_all:bool, drop_dataset:bool, redshift:bool, ):
    """
    """
    if redshift:
        SessionSqlite = base.configure(settings.sqlite_connection_string)
        SessionRedshift = base.configure(settings.redshift_connection_string)

        session = SessionSqlite()
        session_redshift = SessionRedshift()
    else:
        Session = base.configure(settings.sqlite_connection_string)
        session = Session()

    if drop_all:
        logging.info('Dropping the database')
        if redshift:
            session_redshift.execute(redshift_stmts.delete_stmt_all)
            session_redshift.commit()
            Base.metadata.drop_all(bind=session.bind)
        else:
            Base.metadata.drop_all(bind=session.bind)
        
        logging.info('Creating the database')
        if redshift:
            session_redshift.execute(redshift_stmts.create_stmt_measurement)
            session_redshift.commit()
        else:
            Base.metadata.create_all(bind=session.bind)

        gtf = GTFParser(gtf)
        bulk_insert_gtf(session, gtf)

    elif drop_dataset:
        logging.info('Dropping the dataset tables.')
        if redshift:
            session_redshift.execute(redshift_stmts.delete_stmt_measurement)
            session_redshift.commit()
            for n,t in list(base.Base.metadata.tables.items())[::-1]:
                if n in ('sequenceregion','gene','transcript'):
                    continue
                t.drop(bind=session.bind)
        else:
            for n,t in list(base.Base.metadata.tables.items())[::-1]:
                if n in ('sequenceregion','gene','transcript'):
                    continue
                t.drop(bind=session.bind)
        if redshift:
            session.execute(redshift_stmts.create_stmt_measurement)
            session.commit()
        else:
            base.Base.metadata.create_all(bind=session.bind, checkfirst=True)

            

if __name__ == '__main__':
    load_db()