"""test

Revision ID: 3759f4f53660
Revises: 
Create Date: 2024-02-28 16:24:17.085516

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = '3759f4f53660'
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str) -> None:
    globals()["upgrade_%s" % engine_name]()


def downgrade(engine_name: str) -> None:
    globals()["downgrade_%s" % engine_name]()





def upgrade_postgres() -> None:
    pass


def downgrade_postgres() -> None:
    pass


def upgrade_redshift() -> None:
    pass


def downgrade_redshift() -> None:
    pass

