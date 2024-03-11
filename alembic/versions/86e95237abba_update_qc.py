"""update_qc

Revision ID: 86e95237abba
Revises: 9f94af67b82b
Create Date: 2024-03-05 20:19:44.717998

"""
from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = '86e95237abba'
down_revision: Union[str, None] = '9f94af67b82b'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    base.DataSet.set_alembic(revision)
    if engine_name == 'redshift':
        return
    load_db.update_studies_qc(
        connection_string=db_urls['postgres'],
    )

def downgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    if engine_name == 'redshift':
        return
    pass
