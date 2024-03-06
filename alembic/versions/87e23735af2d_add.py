"""add

Revision ID: 87e23735af2d
Revises: 86e95237abba
Create Date: 2024-03-06 19:22:49.228710

"""
from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = '87e23735af2d'
down_revision: Union[str, None] = '86e95237abba'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    if engine_name == 'redshift':
        return
    load_db.add_studies(
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