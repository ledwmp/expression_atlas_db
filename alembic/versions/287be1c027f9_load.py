"""Load

Revision ID: 287be1c027f9
Revises: 9f94af67b82b
Create Date: 2024-03-05 01:28:48.505141

"""
from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = '287be1c027f9'
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    if engine_name == 'redshift':
        return
    load_db.load_db(
        "../veliadb_v0c.gtf",
        drop_all=True,
        use_redshift=True,
        use_s3=True,
        connection_string=db_urls['postgres'],
        redshift_connection_string=db_urls['redshift'],
    )
    load_db.write_studies_qc(
        connection_string=db_urls['postgres'],
    )




def downgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    if engine_name == 'redshift':
        return
    pass
