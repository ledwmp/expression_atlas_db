"""add

Revision ID: 20c0164e17a8
Revises: 3e018f15870d
Create Date: 2024-09-05 23:09:10.386435

"""

from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = "20c0164e17a8"
down_revision: Union[str, None] = "3e018f15870d"
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
