"""update_GSE121212

Revision ID: 163bea0e42b4
Revises: 2fab4bc3bd56
Create Date: 2024-06-05 21:11:47.010961

"""

from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = "163bea0e42b4"
down_revision: Union[str, None] = "2fab4bc3bd56"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    base.DataSet.set_alembic(revision)
    if engine_name == "redshift":
        return
    session = base.configure(db_urls["postgres"])()
    session_redshift = base.configure(db_urls["redshift"])()
    import s3fs

    s3 = s3fs.S3FileSystem()
    study_id = "GSE121212"
    load_db.update_study(
        study_id,
        session,
        session_redshift=session_redshift,
        use_s3=True,
        use_redshift=True,
        s3=s3,
        force=True,
    )
    study = session.query(base.Study).filter(base.Study.internal_id == study_id).all()
    study[0].public = True
    session.commit()
    load_db.write_studies_qc(
        connection_string=db_urls["postgres"],
    )


def downgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    if engine_name == "redshift":
        return
    pass
