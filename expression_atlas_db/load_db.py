import sys
from pathlib import Path
from typing import Dict, List, Union
import logging
import time

import numpy as np
import pandas as pd
import s3fs
from sqlalchemy import Table


def configure_logger(
    log_file=None, level=logging.INFO, overwrite_log=True, format=logging.BASIC_FORMAT
):
    """ """
    # console and file
    if log_file is None:
        logging.basicConfig(stream=sys.stdout, level=level, format=format)
    else:
        logging.basicConfig(
            filename=log_file,
            level=level,
            filemode=("w" if overwrite_log else "a"),
            format=format,
        )
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(format))
        logging.getLogger("").addHandler(console)


configure_logger(
    f"{time.strftime('%Y%m%d_%H%M%S')}_expression_atlas_db.log", level=logging.INFO
)

from expression_atlas_db import base, settings, redshift_stmts, queries
from expression_atlas_db.utils import GTFParser, ExperimentParser, MetaDataFetcher


def bulk_insert_gtf(
    session: base._Session, gtf: Union[GTFParser, None], batch_columns: int = 10000
) -> None:
    """Read from gtf and populate the transcript and gene tables in expression_atlas_db. gtf should
    be an object with transcript_df and gene_df.
    TODO: Replace None above with other sources, perhaps velia_db.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        gtf (Union[GTFParser,None]): object holding transcript_df and gene_df.
        batch_columns (int): number of columns to batch at one time when running inserts.
    """
    # Add all genes from veliadb into sequenceregion, then gene.

    logging.info("Creating gene sequenceregions.")
    records = [
        f"({r.veliadb_id}, '{r.assembly_id}', '{'gene'}')"
        for i, r in gtf.gene_df.iterrows()
    ]
    for i in range(0, int(len(records) / batch_columns) + 1):
        session.execute(
            f"""INSERT INTO sequenceregion (veliadb_id, assembly_id, type) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
        )

    gene_srs = {
        g.veliadb_id: g
        for g in session.query(base.SequenceRegion)
        .filter((base.SequenceRegion.type == "gene"))
        .all()
    }

    logging.info("Populating gene table.")
    records = [
        f"({gene_srs[r.veliadb_id].id}, '{r.gene_id}', '{r.gene_biotype}')"
        for i, r in gtf.gene_df.iterrows()
    ]
    for i in range(0, int(len(records) / batch_columns) + 1):
        session.execute(
            f"""INSERT INTO gene (id, gene_id, gene_biotype) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
        )

    genes = {g.gene_id: g for g in session.query(base.Gene).all()}

    # Add all transcripts from veliadb into sequenceregion, then transcript.

    logging.info("Creating transcript sequence regions.")
    records = [
        f"({r.veliadb_id}, '{r.assembly_id}', '{'transcript'}')"
        for i, r in gtf.transcript_df.iterrows()
    ]
    for i in range(0, int(len(records) / batch_columns) + 1):
        session.execute(
            f"""INSERT INTO sequenceregion (veliadb_id, assembly_id , type) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
        )

    transcript_srs = {
        t.veliadb_id: t
        for t in session.query(base.SequenceRegion)
        .filter((base.SequenceRegion.type == "transcript"))
        .all()
    }

    logging.info("Populating transcript table.")
    records = [
        f"({transcript_srs[r.veliadb_id].id}, '{r.transcript_id}', {genes[r.gene_id].id}, '{r.gene_biotype}')"
        for i, r in gtf.transcript_df.iterrows()
    ]
    for i in range(0, int(len(records) / batch_columns) + 1):
        session.execute(
            f"""INSERT INTO transcript (id, transcript_id, gene_id, gene_biotype) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
        )
    session.commit()


def create_differentialexpression(
    de_df: pd.DataFrame,
    sequenceregions: Dict[str, base.SequenceRegion],
    study: base.Study,
    contrast: base.Contrast,
    batch_columns: int = 10000,
    bulk_insert: bool = False,
    session: Union[base._Session, None] = None,
    de_columns: List[str] = [
        c.name
        for c in base.DifferentialExpression.__table__.columns
        if not c.primary_key and len(c.foreign_keys) == 0
    ],
    s3_fs: Union[s3fs.core.S3FileSystem, None] = None,
    fh: Union[Path, None] = None,
) -> None:
    """Takes a de_df from a contrast and contstructs the table entries via bulk_insert or dumping to a txt file.
    If dumping to text file, either dumps to local if s3_fs is not passed as an argument, or to an s3 bucket.

    Args:
        de_df (pd.DataFrame): DataFrame dumped from ExperimentParser.prepare_differentialexpression.
        sequenceregions (Dict[str,base.SequenceRegion]): Dict of sequenceregions for pulling sequenceregion_ids.
        study (base.Study): Study object.
        contrast (base.Contrast): Contrast object.
        batch_columns (int): Number of columns to batch insert at any time.
        bulk_insert (bool): Run bulk insert or dump table to txt file.
        session (Union[base._Session,None]): SQLAlchemy session object to the main postgres/sqlite db.
        de_coulumns (List[str]): Names of columns from differentialexpresssion table. Used to order txt dumps.
        s3_fs (s3fs.core.S3FileSystem): s3fs object for access to s3.
        fh (Union[Path,None]): Name of file to dump local txt.
    """
    if bulk_insert and not session:
        raise ValueError("Need to provide a session if using bulk insert.")
    if isinstance(s3_fs, s3fs.S3FileSystem) and not fh:
        raise ValueError(
            "Need to pass in s3fs object and s3 location if using s3 dump."
        )
    if not bulk_insert and not fh:
        raise ValueError("Need to pass in fh location if using local dump.")
    if bulk_insert:
        records = [
            f'({contrast.id}, {sequenceregions[i].id}, {", ".join(map(lambda x: str(x),[r[c] for c in de_columns]))})'
            for i, r in de_df.iterrows()
        ]

        logging.info(
            f"Adding de measurements from study {study.velia_id} via bulk insert."
        )
        for i in range(0, int(len(records) / batch_columns) + 1):
            session.execute(
                f"""INSERT INTO differentialexpression (contrast_id, sequenceregion_id, {', '.join(de_columns)}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
            )
            session.commit()
    else:
        records = [
            ",".join(
                map(
                    lambda x: str(x),
                    [contrast.id, sequenceregions[i].id, *[r[c] for c in de_columns]],
                )
            )
            for i, r in de_df.iterrows()
        ]
        logging.info(
            f"Adding de measurements from study {study.velia_id} to csv at {fh}."
        )
        if not isinstance(s3_fs, s3fs.S3FileSystem):
            with open(fh, "w") as f_out:
                f_out.write("\n".join(records))
        else:
            # Check that s3 location exists and is accessible.
            if not s3_fs.exists(str(fh.parent).replace("s3:/", "s3://")):
                raise FileExistsError(
                    f"Location for s3_dump {fh.parent} does not exists."
                )
            with s3_fs.open(str(fh).replace("s3:/", "s3://"), "w") as f_out:
                f_out.write("\n".join(records))
    del records


def create_samplemeasurements(
    m_regions: np.ndarray,
    m_samples: np.ndarray,
    m_measurements: np.ndarray,
    sequenceregions: Dict[str, base.SequenceRegion],
    study: base.Study,
    samples: Dict[str, base.Sample],
    batch_columns: int = 10000,
    bulk_insert: bool = False,
    session: Union[base._Session, None] = None,
    measurement_columns: List[str] = [
        c.name
        for c in base.SampleMeasurement.__table__.columns
        if not c.primary_key and len(c.foreign_keys) == 0
    ],
    s3_fs: Union[s3fs.core.S3FileSystem, None] = None,
    fh: Union[Path, None] = None,
) -> None:
    """Takes the sample measurements from an experiment and create the table entries to samplemeasurement via
    bulk_insert or dumping to a txt file. If dumping to text file, either dumps to local if s3_fs is not passed
    as an argument, or to an s3 bucket.

    Args:
        m_regions (np.ndarray): Array of sequenceregion_ids corresponding to the measurements in m_measurements.
        m_samples (np.ndarray): Array of sample_ids corresponding to the measurements in m_measurements.
        m_measurements (np.ndarray): 2d array of measurements where the last dimension is size as measurement_columns.
        sequenceregions (Dict[str,base.SequenceRegion]): Dict of sequenceregions for pulling sequenceregion_ids.
        study (base.Study): Study object.
        samples (Dict[str,base.Sample]): Dictionary of sample objects.
        batch_columns (int): Number of columns to batch insert at any time.
        bulk_insert (bool): Run bulk insert or dump table to txt file.
        session (Union[base._Session,None]): SQLAlchemy session object to the main postgres/sqlite db.
        measurement_columns (List[str]): Names of columns from differentialexpresssion table. Used to order txt dumps.
        s3_fs (s3fs.core.S3FileSystem): s3fs object for access to s3.
        fh (Union[Path,None]): Name of file to dump local txt.
    """
    if bulk_insert and not session:
        raise ValueError("Need to provide a session if using bulk insert.")
    if isinstance(s3_fs, s3fs.S3FileSystem) and not fh:
        raise ValueError(
            "Need to pass in s3fs object and s3 location if using s3 dump."
        )
    if not bulk_insert and not fh:
        raise ValueError("Need to pass in fh location if using local dump.")
    if bulk_insert:
        records = [
            f'({samples[s].id}, {sequenceregions[i].id}, {", ".join(map(lambda x: str(x), a))})'
            for i, s, a in zip(m_regions, m_samples, m_measurements)
        ]

        logging.info(
            f"Adding sample measurements from study {study.velia_id} via bulk insert."
        )
        for i in range(0, int(len(records) / batch_columns) + 1):
            session.execute(
                f"""INSERT INTO samplemeasurement (sample_id, sequenceregion_id, {', '.join(measurement_columns)}) VALUES 
                {','.join(records[i*batch_columns:(i+1)*batch_columns])};"""
            )
            session.commit()
    else:
        records = [
            ",".join(map(lambda x: str(x), [samples[s].id, sequenceregions[i].id, *a]))
            for i, s, a in zip(m_regions, m_samples, m_measurements)
        ]
        logging.info(
            f"Adding sample measurements from study {study.velia_id} to csv at {fh}."
        )
        if not isinstance(s3_fs, s3fs.S3FileSystem):
            with open(fh, "w") as f_out:
                f_out.write("\n".join(records))
        else:
            # Check that s3 location exists and is accessible.
            if not s3_fs.exists(str(fh.parent).replace("s3:/", "s3://")):
                raise FileExistsError(
                    f"Location for s3_dump {fh.parent} does not exists."
                )
            with s3_fs.open(str(fh).replace("s3:/", "s3://"), "w") as f_out:
                f_out.write("\n".join(records))
    del records


def copy_into_redshift_table(
    session: base._Session,
    table: Table,
    s3_staging_loc: Path,
    iam_role: str = settings.redshift_iam_role,
) -> None:
    """Copies data from s3 bucket into redshift tables using the redshift copy commands.

    Args:
        session (base._Session): SQLAlchemy session object to the redshift db.
        table (base.Base): Table class for destination table.
        s3_staging_loc (Path): The location of the staging directory where txt tables get copied.
        iam_role (str): iam_role string for redshift.
    """
    s3_copy_command = f"""copy {table.name} ({', '.join([c.name for c in table.columns if not c.primary_key])})
    from '{str(s3_staging_loc).replace('s3:/','s3://')}'
    iam_role '{iam_role}'
    delimiter ',';
    """
    try:
        logging.info(f"Copying table {s3_staging_loc} into {table.name}.")
        session.execute(s3_copy_command)
    except Exception as e:
        raise Exception("Unable to copy data into redshift.") from e


def copy_into_db_table(
    session: base._Session,
    table: Table,
    test_staging_loc: Path,
) -> None:
    """Copies data from local into db tables using. If it's ever desired to use a postgres instance
    for samplemeasurment/differentialexpression tables, this needs to be optimized.

    Args:
        session (base._Session): SQLAlchemy session object to the redshift db.
        table (base.Base): Table class for destination table.
        test_staging_loc (Path): The location of the staging directory where txt tables get copied.
    """
    try:
        session.commit()
        logging.info(f"Copying table {test_staging_loc} into {table.name}")
        with open(test_staging_loc, "r") as f_in:
            header = [c.name for c in table.columns if not c.primary_key]
            table_df = pd.read_csv(f_in, header=None, names=header)
        table_df.to_sql(table.name, con=session.bind, index=False, if_exists="append")
    except Exception as e:
        raise Exception("Unable to copy data into db.") from e


def insert_dataset(
    session: base._Session,
    meta: MetaDataFetcher,
    exp: ExperimentParser,
    batch_columns: int = 100000,
    s3: Union[s3fs.core.S3FileSystem, None] = None,
    staging_loc: Path = settings.test_staging_loc,
) -> None:
    """Populates the main db dataset, study, sample, contrast, and samplecontrast tables with the
    data held in meta and exp. Calls create_differentialexpression and create_samplemeasurement to
    dump data to txt files for db upload or run bulk inserts of the data into differentialexpression
    and samplemeasurment tables.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        meta (MetaDataFetcher): Object with metadata for study.
        exp (ExperimentParser): Object with references to study adatas.
        batch_columns (int): Number of columns to batch insert at any time.
        s3 (Union[s3fs.core.S3FileSystem,None]): s3fs object for access to s3.
        staging_loc (Path): Location to stage txt files if not using bulk_insert.
    """

    logging.info(f"Adding study {meta._study_id}.")
    study = base.Study(
        velia_id=meta.velia_id,
        geo_id=meta.geo_id,
        srp_id=meta.srp_id,
        pmid=meta.pmids,
        bio_id=meta.bio_id,
        timestamps=exp.file_timestamps,
        sizes=exp.file_sizes,
        title=meta.project_title,
        description=meta.project_summary,
    )

    session.add(study)

    samples_metadata = meta.samples_metadata
    samples_metadata = samples_metadata.merge(
        exp.samples_metadata[
            exp.samples_metadata.columns[
                ~exp.samples_metadata.columns.isin(samples_metadata.columns)
            ]
        ],
        left_on="Experiment",
        right_index=True,
    ).fillna("NA")

    logging.info(f"Adding {samples_metadata.shape[0]} samples from {meta._study_id}.")
    samples = [
        base.Sample(
            srx_id=s["Experiment"],
            atlas_group=s["sample_type_1"],
            study=study,
            fields=s,
        )
        for s in samples_metadata.to_dict("records")
    ]
    session.add_all(samples)

    sequenceregions = {
        **{g.gene_id: g for g in session.query(base.Gene).all()},
        **{t.transcript_id: t for t in session.query(base.Transcript).all()},
    }

    measurement_columns = [
        base.SampleMeasurement._column_map[c.name]
        for c in base.SampleMeasurement.__table__.columns
        if not c.primary_key and len(c.foreign_keys) == 0
    ]

    de_columns = [
        base.DifferentialExpression._column_map[c.name]
        for c in base.DifferentialExpression.__table__.columns
        if not c.primary_key and len(c.foreign_keys) == 0
    ]

    gm_samples, gm_regions, gm_measurments, g_measurements = exp.prepare_measurements(
        measurement_type="gene", measurements=measurement_columns
    )
    g_exp_dict = exp.prepare_differentialexpression(
        measurement_type="gene", de_columns=de_columns
    )
    tm_samples, tm_regions, tm_measurments, t_measurements = exp.prepare_measurements(
        measurement_type="transcript", measurements=measurement_columns
    )
    t_exp_dict = exp.prepare_differentialexpression(
        measurement_type="transcript", de_columns=de_columns
    )

    for c, (c_df, c_left, c_right) in g_exp_dict.items():

        logging.info(f"Adding contrast {c} from study {study.velia_id}.")

        left_condition_display = c_df["left_condition_display"].unique()[0]
        right_condition_display = c_df["right_condition_display"].unique()[0]
        sample_condition_key = c_df["sample_condition_key"].unique()[0]
        left_condition_level = c_df["left_condition_level"].unique()[0]
        right_condition_level = c_df["right_condition_level"].unique()[0]

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

        logging.info(
            f"Adding left samplecontrasts from contrast {c} from study {study.velia_id}."
        )
        left_samples = session.query(base.Sample).filter(base.Sample.srx_id.in_(c_left))
        left_sample_contrasts = [
            base.SampleContrast(
                sample=s,
                contrast_side="left",
                contrast=contrast,
                study=study,
                condition_display=contrast.left_condition_display,
                condition_level=contrast.left_condition_level,
                sample_condition_key=contrast.sample_condition_key,
            )
            for s in left_samples
        ]

        logging.info(
            f"Adding right samplecontrasts from contrast {c} from study {study.velia_id}."
        )
        right_samples = session.query(base.Sample).filter(
            base.Sample.srx_id.in_(c_right)
        )
        right_sample_contrasts = [
            base.SampleContrast(
                sample=s,
                contrast_side="right",
                contrast=contrast,
                study=study,
                condition_display=contrast.right_condition_display,
                condition_level=contrast.right_condition_level,
                sample_condition_key=contrast.sample_condition_key,
            )
            for s in right_samples
        ]
        session.add_all(left_sample_contrasts + right_sample_contrasts)

        create_differentialexpression(
            c_df,
            sequenceregions,
            study,
            contrast,
            s3_fs=s3,
            fh=staging_loc / f"gene_de.{study.id}.{contrast.id}.csv",
        )

    samples = {
        s.srx_id: s
        for s in session.query(base.Sample).filter(
            base.Sample.srx_id.in_(np.unique(gm_samples))
        )
    }
    create_samplemeasurements(
        gm_regions,
        gm_samples,
        gm_measurments,
        sequenceregions,
        study,
        samples,
        s3_fs=s3,
        fh=staging_loc / f"gene_measurement.{study.id}.csv",
    )

    for c, (c_df, c_left, c_right) in t_exp_dict.items():

        contrast = (
            session.query(base.Contrast)
            .filter(
                (base.Contrast.contrast_name == c)
                & (base.Contrast.study_id == study.id)
            )
            .first()
        )

        create_differentialexpression(
            c_df,
            sequenceregions,
            study,
            contrast,
            s3_fs=s3,
            fh=staging_loc / f"transcript_de.{study.id}.{contrast.id}.csv",
        )

    samples = {
        s.srx_id: s
        for s in session.query(base.Sample).filter(
            base.Sample.srx_id.in_(np.unique(tm_samples))
        )
    }
    create_samplemeasurements(
        tm_regions,
        tm_samples,
        tm_measurments,
        sequenceregions,
        study,
        samples,
        s3_fs=s3,
        fh=staging_loc / f"transcript_measurement.{study.id}.csv",
    )
    del gm_samples
    del gm_regions
    del gm_measurments
    del tm_samples
    del tm_regions
    del tm_measurments


def delete_studies_from_update(
    alembic_revision_id: str,
    use_redshift: bool = False,
    use_s3: bool = False,
    connection_string: str = settings.db_connection_string,
    redshift_connection_string: Union[str, None] = settings.redshift_connection_string,
) -> None:
    """Removes all data associated with a specific alembic update from db. If data were inserted into
    redshift tables via copy, deletes staged files in s3 located in staging loc. Deletes all studies from
    db with specifc alembic id which cascade deletes all entries in dataset table.

    Args:
        alembic_revision_id (str): alembic_revision_id to be deleted from db.
        use_redshift (bool): Delete out of redshift samplemeasurement and differentialexpression tables.
        use_s3 (bool): Delete staged files out of s3 bucket.
        connection_string (str): Connection string to main postgres/sqlite db.
        redshift_connection_string (str): Connection string to redshift db.
    """
    logging.info(
        f"Querying studies across atlas with revision_id: {alembic_revision_id}."
    )

    if use_redshift:
        Session = base.configure(connection_string)
        SessionRedshift = base.configure(redshift_connection_string)

        session = Session()
        session_redshift = SessionRedshift()

    else:
        Session = base.configure(connection_string)
        session = Session()

    s3 = s3fs.S3FileSystem() if use_s3 else None

    velia_ids = [
        e
        for (e,) in session.query(base.Study.velia_id).filter(
            base.Study.alembic_id == alembic_revision_id
        )
    ]

    logging.info(
        f"Found {len(velia_ids)} studies to be deleted with revision_id: {alembic_revision_id}."
    )

    for e in set(velia_ids):

        try:
            delete_study(
                e,
                session,
                session_redshift=session_redshift if use_redshift else None,
                use_s3=use_s3,
                use_redshift=use_redshift,
                s3=s3 if use_s3 else None,
            )
        except Exception as e:
            logging.exception(e)
            session.rollback()
            if use_redshift:
                session_redshift.rollback()
    session.close()
    if use_redshift:
        session_redshift.close()


def delete_study(
    velia_study: str,
    session: base._Session,
    session_redshift: Union[base._Session, None] = None,
    use_s3: bool = False,
    use_redshift: bool = False,
    s3: Union[s3fs.core.S3FileSystem, None] = None,
) -> None:
    """Deletes a specific velia_study out of db. If data were inserted into redshift tables via
    copy, deletes staged files in s3 located in staging loc.

    Args:
        velia_study (str): velia_study id, usually a geo or srp id.
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        use_s3 (bool): Delete staged files out of s3 bucket.
        use_redshift (bool): Delete out of redshift samplemeasurement and differentialexpression tables.
        s3 (s3fs.core.s3FileSysetem): s3fs object for access to s3.
    """
    logging.info(f"Deleting study: {velia_study} from expression_atlas_db.")

    study = session.query(base.Study).filter(base.Study.velia_id == velia_study).all()

    # Delete the sample_measurements/differential_expression files out of s3_staging_loc.
    gene_measurement_loc = str(
        Path(settings.s3_staging_loc if use_s3 else settings.test_staging_loc)
        / f"gene_measurement.{study[0].id}.csv"
    ).replace("s3:/", "s3://")
    transcript_measurement_loc = str(
        Path(settings.s3_staging_loc if use_s3 else settings.test_staging_loc)
        / f"transcript_measurement.{study[0].id}.csv"
    ).replace("s3:/", "s3://")

    logging.info(f"Removing staged file at: {gene_measurement_loc}.")
    if use_s3:
        s3.rm(gene_measurement_loc)
    else:
        Path(gene_measurement_loc).unlink()
    logging.info(f"Removing staged file at: {transcript_measurement_loc}.")
    if use_s3:
        s3.rm(transcript_measurement_loc)
    else:
        Path(transcript_measurement_loc).unlink()

    contrasts = (
        session.query(base.Contrast).filter(base.Contrast.study_id == study[0].id).all()
    )

    for c in contrasts:
        gene_de_loc = str(
            Path(settings.s3_staging_loc if use_s3 else settings.test_staging_loc)
            / f"gene_de.{study[0].id}.{c.id}.csv"
        ).replace("s3:/", "s3://")
        transcript_de_loc = str(
            Path(settings.s3_staging_loc if use_s3 else settings.test_staging_loc)
            / f"transcript_de.{study[0].id}.{c.id}.csv"
        ).replace("s3:/", "s3://")
        logging.info(f"Removing staged file at: {gene_de_loc}.")
        if use_s3:
            s3.rm(gene_de_loc)
        else:
            Path(gene_de_loc).unlink()
        logging.info(f"Removing staged file at: {transcript_de_loc}.")
        if use_s3:
            s3.rm(transcript_de_loc)
        else:
            Path(transcript_de_loc).unlink()

    samples = (
        session.query(base.Sample).filter(base.Sample.study_id == study[0].id).all()
    )

    # Delete the sample_measurements/differential_expression entries out of redshift tables.
    if use_redshift:
        logging.info(f"Removing entries out of redshift differentialexpression table.")
        session_redshift.execute(
            f"DELETE FROM differentialexpression WHERE contrast_id in ({', '.join([str(c.id) for c in contrasts])});"
        )
        logging.info(f"Removing entries out of redshift samplemeasurement table.")
        session_redshift.execute(
            f"DELETE FROM samplemeasurement WHERE sample_id in ({', '.join([str(s.id) for s in samples])});"
        )
    else:
        logging.info(
            f"Removing entries out of local/postgres differentialexpression table."
        )
        session.execute(
            f"DELETE FROM differentialexpression WHERE contrast_id in ({', '.join([str(c.id) for c in contrasts])});"
        )
        logging.info(f"Removing entries out of local/postgres samplemeasurement table.")
        session.execute(
            f"DELETE FROM samplemeasurement WHERE sample_id in ({', '.join([str(s.id) for s in samples])});"
        )

    # Delete the study, should cascade and delete all samples, contrasts, and samplecontrasts from dataset.
    logging.info(f"Removing study: {velia_study} from expression_atlas_db.")
    session.delete(study[0])
    session.commit()
    if use_redshift:
        session_redshift.commit()
    logging.info(f"Fully deleted study: {velia_study}.")


def update_study(
    velia_study: str,
    session: base._Session,
    session_redshift: Union[base._Session, None] = None,
    use_s3: bool = False,
    use_redshift: bool = False,
    update_timestamp: bool = False,
    s3: Union[s3fs.core.S3FileSystem, None] = None,
    force: bool = False,
) -> None:
    """Updates a specific velia_study present in db. If data were inserted into redshift tables via
    copy, deletes staged files in s3 located in staging loc and replaces with new staged files. Stats
    adatas locally or in s3 to check for changes since last update. If file sizes have changed, deletes
    studies from db and re-inserts updated data. If update_timestamp flag is set, also updates upon changes to
    the file timestamp.

    Args:
        velia_study (str): velia_study id, usually a geo or srp id.
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        use_s3 (bool): Delete and replaced staged files in s3 bucket.
        use_redshift (bool): Delete and replace rows in redshift samplemeasurement and differentialexpression tables.
        update_timestamp (bool): Default is to just assess files to update based on file sizes, this also checks timestamps.
        s3 (s3fs.core.s3FileSysetem): s3fs object for access to s3.
        force (bool) Force update.
    """
    logging.info(f"Checking study for update: {velia_study}.")
    exp = ExperimentParser(
        velia_study,
        Path(settings.s3_experiment_loc if use_s3 else settings.test_experiment_loc),
    )
    if use_s3:
        exp.enable_s3(s3)

    # Check the size and date of access against the size recorded in velia_db.
    exp.stat_adatas()

    study = session.query(base.Study).filter(base.Study.velia_id == velia_study).all()

    if len(study) != 1:
        raise ValueError(
            "Number of studies to update > 1. Study duplicated or does not exist in db."
        )

    if study[0].sizes != exp.file_sizes:
        logging.info(
            f"Updating: {velia_study} file_sizes mismatch {study[0].sizes} -> {exp.file_sizes}."
        )
    elif update_timestamp and study[0].timestamps != exp.file_timestamps:
        logging.info(
            f"Updating: {velia_study} timestamps mismatch {study[0].timestamps} -> {exp.file_timestamps}."
        )
    elif force:
        logging.info(f"Updating: {velia_study} forcefully.")
    else:
        logging.info(f"Nothing to update: {velia_study}")
        return

    delete_study(
        velia_study,
        session,
        session_redshift,
        use_s3=use_s3,
        use_redshift=use_redshift,
        s3=s3,
    )

    logging.info(f"Re-loading adatas: {velia_study}.")
    exp.load_adatas()
    logging.info(f"Re-fetching metadata: {velia_study}.")
    meta = MetaDataFetcher(velia_study, exp.samples)
    logging.info(f"Re-inserting datasets: {velia_study}.")
    insert_dataset(
        session,
        meta,
        exp,
        s3=s3,
        staging_loc=Path(
            settings.s3_staging_loc if use_s3 else settings.test_staging_loc
        ),
    )
    study_id = [
        *session.query(base.Study.id).filter(base.Study.velia_id == velia_study)
    ][0]
    contrast_ids = [
        *session.query(base.Contrast.id).filter(
            base.Contrast.study.has(velia_id=velia_study)
        )
    ]
    if use_redshift:
        for m in ("gene", "transcript"):
            copy_into_redshift_table(
                session_redshift,
                base.Base.metadata.tables["samplemeasurement"],
                Path(settings.s3_staging_loc) / f"{m}_measurement.{study_id[0]}.csv",
            )
            for c in contrast_ids:
                copy_into_redshift_table(
                    session_redshift,
                    base.Base.metadata.tables["differentialexpression"],
                    Path(settings.s3_staging_loc) / f"{m}_de.{study_id[0]}.{c[0]}.csv",
                )
    else:
        for m in ("gene", "transcript"):
            copy_into_db_table(
                session,
                base.Base.metadata.tables["samplemeasurement"],
                Path(settings.test_staging_loc) / f"{m}_measurement.{study_id[0]}.csv",
            )
            for c in contrast_ids:
                copy_into_db_table(
                    session,
                    base.Base.metadata.tables["differentialexpression"],
                    Path(settings.test_staging_loc)
                    / f"{m}_de.{study_id[0]}.{c[0]}.csv",
                )
    # Will need to find a way to catch exceptions in one database and rollback in the other.
    session.commit()
    if use_redshift:
        session_redshift.commit()
    del exp
    del meta


def update_studies(
    use_redshift: bool = False,
    use_s3: bool = False,
    connection_string: str = settings.db_connection_string,
    redshift_connection_string: Union[str, None] = settings.redshift_connection_string,
    update_timestamp: bool = False,
) -> None:
    """Checks for files to update in staging_loc based upon data in db. Stats adatas locally or in s3 to check
    for changes since last update. If file sizes have changed, deletes studies from db and re-inserts updated data.
    If update_timestamp flag is set, also updates upon changes to the file timestamp.

    Args:
        use_redshift (bool): Delete and replace rows in redshift samplemeasurement and differentialexpression tables.
        use_s3 (bool): Delete and replaced staged files in s3 bucket.
        connection_string (str): Connection string to main postgres/sqlite db.
        redshift_connection_string (str): Connection string to redshift db.
        update_timestamp (bool): Default is to just assess files to update based on file sizes, this also checks timestamps.
    """
    if use_redshift:
        Session = base.configure(connection_string)
        SessionRedshift = base.configure(redshift_connection_string)

        session = Session()
        session_redshift = SessionRedshift()

    else:
        Session = base.configure(connection_string)
        session = Session()

    s3 = s3fs.S3FileSystem() if use_s3 else None
    if use_s3:
        velia_ids = [
            Path(f).parts[-2]
            for f in s3.glob(
                str(Path(settings.s3_experiment_loc, "./*/de_results/")).replace(
                    "s3:/", "s3://"
                )
            )
        ]
    else:
        velia_ids = [
            p for p in Path(settings.test_experiment_loc).iterdir() if p.is_dir()
        ]

    logging.info(
        f"Found {len(velia_ids)} at {settings.s3_experiment_loc if use_s3 else settings.test_experiment_loc}."
    )

    existing_studies = [e for (e,) in session.query(base.Study.velia_id)]

    logging.info(f"Found {len(existing_studies)} already populated in database.")

    for e in set(velia_ids).difference(existing_studies):
        logging.info(f"Skipping: {e} does not exist in database.")

    for e in set(velia_ids).intersection(existing_studies):

        try:
            update_study(
                e,
                session,
                session_redshift=session_redshift if use_redshift else None,
                use_s3=use_s3,
                use_redshift=use_redshift,
                s3=s3 if use_s3 else None,
                update_timestamp=update_timestamp,
            )
        except Exception as e:
            logging.exception(e)
            session.rollback()
            if use_redshift:
                session_redshift.rollback()
    session.close()
    if use_redshift:
        session_redshift.close()


def update_studies_qc(
    connection_string: str = settings.db_connection_string,
    qc_loc: Path = Path(settings.s3_staging_loc),
) -> None:
    """Pulls qc file from s3 bucket and updates public and quality columns in db.

    Args:
        connection_string (str): Connection string to main postgres/sqlite db.
        qc_loc (Path): Location in s3 for qc files (defaults to the staging_loc).
    """
    Session = base.configure(connection_string)
    session = Session()

    logging.info("Updating QC sheet.")

    s3 = s3fs.S3FileSystem()

    qc_files = s3.glob(str(qc_loc / "qc*txt").replace("s3:/", "s3://"))
    qc_files = sorted(
        qc_files, key=lambda x: str(Path(x).parts[-1].split(".")[1]), reverse=True
    )

    logging.info(f"Reading QC sheet: {qc_files[0]}.")

    with s3.open(qc_files[0], "rb") as f_in:
        qc_df = pd.read_csv(f_in, sep="|")

    existing_studies = [s for (s,) in session.query(base.Study.velia_id)]

    if len(set(qc_df["velia_id"]).difference(existing_studies)) > 0:
        raise ValueError("Study found in QC document not populated in database.")

    for r in qc_df.to_dict("records"):
        session.query(base.Study).filter(base.Study.velia_id == r["velia_id"]).update(r)

    session.commit()
    logging.info(f"Updated studies with QC sheet: {qc_files[0]}.")


def write_studies_qc(
    connection_string: str = settings.db_connection_string,
    qc_loc: Path = Path(settings.s3_staging_loc),
) -> None:
    """Dumps a new qc template to qc_loc with current studies in db and their public/quality fields.

    Args:
        connection_string (str): Connection string to main postgres/sqlite db.
        qc_loc (Path): Location in s3 for qc files (defaults to the staging_loc).
    """
    Session = base.configure(connection_string)
    session = Session()

    logging.info("Writing QC sheet.")

    s3 = s3fs.S3FileSystem()

    qc_files = s3.glob(str(qc_loc / "qc*txt").replace("s3:/", "s3://"))
    qc_files = sorted(
        qc_files, key=lambda x: str(Path(x).parts[-1].split(".")[1]), reverse=True
    )
    if len(qc_files) > 0:
        qc_number = int(Path(qc_files[0]).parts[-1].split(".")[1]) + 1
        logging.info(
            f"Found last QC sheet: {qc_files[0]}, updating to number: {qc_number}."
        )
    else:
        qc_number = 0
        logging.info(f"No QC sheets found, using number: {qc_number}.")

    with s3.open(
        str(Path(qc_loc) / f"qc.{qc_number}.txt").replace("s3:/", "s3://"), "w"
    ) as f_out:
        queries.fetch_studies(session, public=False)[
            ["velia_id", "srp_id", "geo_id", "public", "quality"]
        ].to_csv(f_out, index=False, sep="|")

    logging.info("QC sheet updated.")


def add_study(
    velia_study: str,
    session: base._Session,
    session_redshift: Union[base._Session, None] = None,
    use_s3: bool = False,
    use_redshift: bool = False,
    s3: Union[s3fs.core.S3FileSystem, None] = None,
) -> None:
    """Adds a specific velia_study to db. Parses experiment and fetches metadata related to experiment.
    Calls insert_dataset to add studies to postgres/sqlite database. If use_s3, dumps files into staging_loc
    that get uploaded into redshift db, otherwise dumps local txt files or runs bulk_insert.

    Args:
        velia_study (str): velia_study id, usually a geo or srp id.
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        use_s3 (bool): Delete staged files out of s3 bucket.
        use_redshift (bool): Delete out of redshift samplemeasurement and differentialexpression tables.
        s3 (s3fs.core.s3FileSysetem): s3fs object for access to s3.
    """

    logging.info(f"Parsing: {velia_study}.")
    exp = ExperimentParser(
        velia_study,
        Path(settings.s3_experiment_loc if use_s3 else settings.test_experiment_loc),
    )
    if use_s3:
        exp.enable_s3(s3)
    exp.load_adatas()
    logging.info(f"Fetching metadata: {velia_study}.")
    meta = MetaDataFetcher(velia_study, exp.samples)
    logging.info(f"Inserting datasets: {velia_study}.")
    insert_dataset(
        session,
        meta,
        exp,
        s3=s3,
        staging_loc=Path(
            settings.s3_staging_loc if use_s3 else settings.test_staging_loc
        ),
    )
    study_id = [
        *session.query(base.Study.id).filter(base.Study.velia_id == velia_study)
    ][0]
    contrast_ids = [
        *session.query(base.Contrast.id).filter(
            base.Contrast.study.has(velia_id=velia_study)
        )
    ]
    if use_redshift:
        for m in ("gene", "transcript"):
            copy_into_redshift_table(
                session_redshift,
                base.Base.metadata.tables["samplemeasurement"],
                Path(settings.s3_staging_loc) / f"{m}_measurement.{study_id[0]}.csv",
            )
            for c in contrast_ids:
                copy_into_redshift_table(
                    session_redshift,
                    base.Base.metadata.tables["differentialexpression"],
                    Path(settings.s3_staging_loc) / f"{m}_de.{study_id[0]}.{c[0]}.csv",
                )
    else:
        for m in ("gene", "transcript"):
            copy_into_db_table(
                session,
                base.Base.metadata.tables["samplemeasurement"],
                Path(settings.test_staging_loc) / f"{m}_measurement.{study_id[0]}.csv",
            )
            for c in contrast_ids:
                copy_into_db_table(
                    session,
                    base.Base.metadata.tables["differentialexpression"],
                    Path(settings.test_staging_loc)
                    / f"{m}_de.{study_id[0]}.{c[0]}.csv",
                )

    # Will need to find a way to catch exceptions in one database and rollback in the other.
    session.commit()
    if use_redshift:
        session_redshift.commit()
    del exp
    del meta


def add_studies(
    use_redshift: bool = False,
    use_s3: bool = False,
    connection_string: str = settings.db_connection_string,
    redshift_connection_string: Union[str, None] = settings.redshift_connection_string,
) -> None:
    """Adds studies located in experiment_loc to db. Checks db for existence of study, then calls add_study
    to populate study into db.

    Args:
        use_redshift (bool): Add to redshift samplemeasurement and differentialexpression tables.
        use_s3 (bool): Stage files in s3 bucket, pull experiments from s3_experiment_loc.
        connection_string (str): Connection string to main postgres/sqlite db.
        redshift_connection_string (str): Connection string to redshift db.
    """
    if use_redshift:
        Session = base.configure(connection_string)
        SessionRedshift = base.configure(redshift_connection_string)

        session = Session()
        session_redshift = SessionRedshift()

    else:
        Session = base.configure(connection_string)
        session = Session()

    s3 = s3fs.S3FileSystem() if use_s3 else None
    if use_s3:
        velia_ids = [
            Path(f).parts[-2]
            for f in s3.glob(
                str(Path(settings.s3_experiment_loc, "./*/de_results/")).replace(
                    "s3:/", "s3://"
                )
            )
        ]
    else:
        velia_ids = [
            p for p in Path(settings.test_experiment_loc).iterdir() if p.is_dir()
        ]

    logging.info(
        f"Found {len(velia_ids)} at {settings.s3_experiment_loc if use_s3 else settings.test_experiment_loc}."
    )

    existing_studies = [e for (e,) in session.query(base.Study.velia_id)]

    logging.info(f"Found {len(existing_studies)} already populated in database.")

    for e in set(velia_ids).intersection(existing_studies):
        logging.info(f"Skipping: {e} already exists in database.")

    for e in set(velia_ids).difference(existing_studies):

        try:
            add_study(
                e,
                session,
                session_redshift=session_redshift if use_redshift else None,
                use_s3=use_s3,
                use_redshift=use_redshift,
                s3=s3 if use_s3 else None,
            )
        except Exception as e:
            logging.exception(e)
            session.rollback()
            if use_redshift:
                session_redshift.rollback()
    session.close()
    if use_redshift:
        session_redshift.close()


def load_db(
    gtf: str,
    drop_all: bool = False,
    drop_dataset: bool = False,
    use_redshift: bool = False,
    use_s3: bool = False,
    connection_string: str = settings.db_dev_connection_string,
    redshift_connection_string: Union[
        str, None
    ] = settings.redshift_dev_connection_string,
) -> None:
    """Initial load_db method. If drop_all, calls bulk_insert_gtf to populate tables inheriting from
    sequenceregion table. Calls add_study on all experiments located at experiment_loc and populates
    dataset tables.

    Args:
        gtf (str): Path to gtf to load into sequenceregions, gene, and transcript tables.
        drop_all (bool): Drop all tables.
        drop_dataset (bool): Keep sequenceregions tables, just drop tables inheriting from dataset.
        use_redshift (bool): Add to redshift samplemeasurement and differentialexpression tables.
        use_s3 (bool): Stage files in s3 bucket, pull experiments from s3_experiment_loc.
        connection_string (str): Connection string to main postgres/sqlite db.
        redshift_connection_string (str): Connection string to redshift db.
    """
    if use_redshift:
        Session = base.configure(connection_string)
        SessionRedshift = base.configure(redshift_connection_string)

        session = Session()
        session_redshift = SessionRedshift()
    else:
        Session = base.configure(connection_string)
        session = Session()

    if drop_all:
        logging.info("Dropping the database.")
        if use_redshift:
            session_redshift.execute(redshift_stmts.delete_stmt_all)
            session_redshift.commit()
            base.Base.metadata.drop_all(bind=session.bind)
        else:
            base.Base.metadata.drop_all(bind=session.bind)

        logging.info("Creating the database.")
        if use_redshift:
            session_redshift.execute(redshift_stmts.create_stmt_measurement)
            session_redshift.commit()
            base.Base.metadata.create_all(bind=session.bind)
        else:
            base.Base.metadata.create_all(bind=session.bind)

        gtf = GTFParser(gtf)
        bulk_insert_gtf(session, gtf)
        session.commit()

    elif drop_dataset:
        logging.info("Dropping the dataset tables.")
        if use_redshift:
            session_redshift.execute(redshift_stmts.delete_stmt_all)
            session_redshift.commit()
            for n, t in list(base.Base.metadata.tables.items())[::-1]:
                if n in ("sequenceregion", "gene", "transcript"):
                    continue
                t.drop(bind=session.bind, checkfirst=True)
        else:
            for n, t in list(base.Base.metadata.tables.items())[::-1]:
                if n in ("sequenceregion", "gene", "transcript"):
                    continue
                t.drop(bind=session.bind, checkfirst=True)
        if use_redshift:
            session_redshift.execute(redshift_stmts.create_stmt_measurement)
            session_redshift.commit()
            base.Base.metadata.create_all(bind=session.bind, checkfirst=True)
        else:
            base.Base.metadata.create_all(bind=session.bind, checkfirst=True)

    s3 = s3fs.S3FileSystem() if use_s3 else None
    if use_s3:
        velia_studies = [
            Path(f).parts[-2]
            for f in s3.glob(
                str(Path(settings.s3_experiment_loc, "./*/de_results/")).replace(
                    "s3:/", "s3://"
                )
            )
        ]
    else:
        velia_studies = [
            p.stem for p in Path(settings.test_experiment_loc).iterdir() if p.is_dir()
        ]

    logging.info(f"Loading {len(velia_studies)} studies into expression_atlas_db.")
    for e in velia_studies:
        try:
            add_study(
                e,
                session,
                session_redshift=session_redshift if use_redshift else None,
                use_s3=use_s3,
                use_redshift=use_redshift,
                s3=s3 if use_s3 else None,
            )
        except Exception as e:
            logging.exception(e)
            session.rollback()
            if use_redshift:
                session_redshift.rollback()
    session.close()
    if use_redshift:
        session_redshift.close()


if __name__ == "__main__":
    pass
