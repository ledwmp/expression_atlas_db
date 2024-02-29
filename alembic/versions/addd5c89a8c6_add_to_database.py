"""Add to database

Revision ID: addd5c89a8c6
Revises: 
Create Date: 2024-02-29 21:30:48.076543

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db

# revision identifiers, used by Alembic.
revision: str = 'addd5c89a8c6'
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str) -> None:
    if engine_name == 'redshift':
        return
    load_db.add_studies(use_redshift=True, use_s3=True)


def downgrade(engine_name: str) -> None:
    if engine_name == 'redshift':
        return
    pass
