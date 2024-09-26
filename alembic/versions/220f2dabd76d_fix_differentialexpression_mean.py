"""fix_differentialexpression_mean

Revision ID: 220f2dabd76d
Revises: 4ce112652c16
Create Date: 2024-09-26 01:11:12.190841

"""

from pathlib import Path
from typing import Sequence, Union, Dict
import logging

from alembic import op
import sqlalchemy as sa
import s3fs

from expression_atlas_db import base, load_db, settings, utils

# revision identifiers, used by Alembic.
revision: str = "220f2dabd76d"
down_revision: Union[str, None] = "4ce112652c16"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None

stmt_create_staging_table = """CREATE TEMP TABLE stage (
	id INTEGER NOT NULL IDENTITY(1,1), 
	contrast_id INTEGER, 
	sequenceregion_id INTEGER, 
	basemean FLOAT, 
	log2foldchange FLOAT, 
	lfcse FLOAT, 
	stat FLOAT, 
	pvalue DOUBLE PRECISION, 
	log10_pvalue FLOAT,
	padj DOUBLE PRECISION, 
	log10_padj FLOAT, 
	control_mean FLOAT, 
	case_mean FLOAT, 
	PRIMARY KEY (id)
);"""

s3_copy_command = f"""COPY stage ({', '.join([c.name for c in base.DifferentialExpression.__table__.columns if not c.primary_key])})
FROM '{{table_fh}}'
iam_role '{settings.redshift_iam_role}'
delimiter ',';"""

stmt_merge_update = """BEGIN;
LOCK differentialexpression;
UPDATE differentialexpression SET 
case_mean = stage.case_mean, 
control_mean = stage.control_mean 
FROM stage 
WHERE differentialexpression.contrast_id = stage.contrast_id 
and differentialexpression.sequenceregion_id = stage.sequenceregion_id;
DROP TABLE stage;
END;"""


def update_differential_expression_table(
    session: base._Session,
    session_redshift: base._Session,
    velia_id: str,
    s3: Union[s3fs.S3FileSystem, None] = s3fs.S3FileSystem(),
    staging_loc: Path = Path(settings.s3_staging_loc),
) -> None:
    """ """
    study = session.query(base.Study).filter(base.Study.velia_id == velia_id).all()
    if len(study) != 1:
        raise Exception(f"More than one or no study found for velia_id: {velia_id}.")
    study = study[0]

    sequenceregions = {
        **{g.gene_id: g for g in session.query(base.Gene).all()},
        **{t.transcript_id: t for t in session.query(base.Transcript).all()},
    }

    logging.info(f"Reading adatas {velia_id}...")
    exp = utils.ExperimentParser(velia_id, Path(settings.s3_experiment_loc))
    exp.enable_s3(s3)
    exp.load_adatas()

    # Create column names for differentialexpression table.

    de_columns = [
        base.DifferentialExpression._column_map[c.name]
        for c in base.DifferentialExpression.__table__.columns
        if not c.primary_key and len(c.foreign_keys) == 0
    ]

    # Load gene contrasts.

    g_exp_dict = exp.prepare_differentialexpression(
        measurement_type="gene",
        de_columns=de_columns,
    )

    # Load transcript contrasts.

    t_exp_dict = exp.prepare_differentialexpression(
        measurement_type="transcript",
        de_columns=de_columns,
    )

    # Iterate through gene contrasts. Push upadated table, correct entries in differentialexpression table.

    for c, (c_df, _, _) in g_exp_dict.items():

        logging.info(f"Modifying contrast {c} from study {study.velia_id}.")

        contrast = (
            session.query(base.Contrast)
            .filter(base.Contrast.contrast_name == c)
            .filter(base.Contrast.study_id == study.id)
            .first()
        )

        table_fh = staging_loc / f"gene_de.{study.id}.{contrast.id}.csv"

        load_db.create_differentialexpression(
            c_df,
            sequenceregions,
            study,
            contrast,
            fh=table_fh,
            s3_fs=s3,
        )

        # Load/replace.

        session_redshift.execute(stmt_create_staging_table)
        session_redshift.execute(
            s3_copy_command.format(table_fh=str(table_fh).replace("s3:/", "s3://"))
        )
        session_redshift.execute(stmt_merge_update)

    # Iterate through transcript contrasts. Push upadated table, correct entries in differentialexpression table.

    for c, (c_df, _, _) in t_exp_dict.items():

        logging.info(f"Modifying contrast {c} from study {study.velia_id}.")

        contrast = (
            session.query(base.Contrast)
            .filter(base.Contrast.contrast_name == c)
            .filter(base.Contrast.study_id == study.id)
            .first()
        )

        table_fh = staging_loc / f"transcript_de.{study.id}.{contrast.id}.csv"

        load_db.create_differentialexpression(
            c_df,
            sequenceregions,
            study,
            contrast,
            fh=table_fh,
            s3_fs=s3,
        )

        # Load/replace.

        session_redshift.execute(stmt_create_staging_table)
        session_redshift.execute(
            s3_copy_command.format(table_fh=str(table_fh).replace("s3:/", "s3://"))
        )
        session_redshift.execute(stmt_merge_update)


def upgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    base.DataSet.set_alembic(revision)
    if engine_name == "redshift":
        return

    s3 = s3fs.S3FileSystem()
    velia_ids = [
        Path(f).parts[-2]
        for f in s3.glob(
            str(Path(settings.s3_experiment_loc, "./*/de_results/")).replace(
                "s3:/", "s3://"
            )
        )
    ]

    session = base.configure(db_urls["postgres"])()
    session_redshift = base.configure(db_urls["redshift"])()

    for vid in velia_ids:
        try:
            update_differential_expression_table(
                session,
                session_redshift,
                vid,
                s3=s3,
            )
        except Exception as e:
            logging.exception(e)
            session.rollback()
            session_redshift.rollback()


def downgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    if engine_name == "redshift":
        return
    pass
