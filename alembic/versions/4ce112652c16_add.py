"""add

Revision ID: 4ce112652c16
Revises: cdf22046d53d
Create Date: 2024-09-12 15:06:09.616696

"""

from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = "4ce112652c16"
down_revision: Union[str, None] = "cdf22046d53d"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    base.DataSet.set_alembic(revision)
    if engine_name == "redshift":
        return
    load_db.add_studies(
        use_redshift=True,
        use_s3=True,
        connection_string=db_urls["postgres"],
        redshift_connection_string=db_urls["redshift"],
    )
    load_db.write_studies_qc(
        connection_string=db_urls["postgres"],
    )


def downgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    if engine_name == "redshift":
        return
    pass