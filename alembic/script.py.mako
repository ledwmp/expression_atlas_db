<%!
import re

%>"""${message}

Revision ID: ${up_revision}
Revises: ${down_revision | comma,n}
Create Date: ${create_date}

"""
from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa
${imports if imports else ""}

from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = ${repr(up_revision)}
down_revision: Union[str, None] = ${repr(down_revision)}
branch_labels: Union[str, Sequence[str], None] = ${repr(branch_labels)}
depends_on: Union[str, Sequence[str], None] = ${repr(depends_on)}


def upgrade(engine_name: str, db_urls: Dict[str,str]) -> None:
    base.DataSet.set_alembic(revision)
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
