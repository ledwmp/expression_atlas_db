"""add_studyqueue_table

Revision ID: c6f2acee9b09
Revises: 3726b51b1db8
Create Date: 2024-07-29 16:31:26.992971

"""

from typing import Sequence, Union, Dict

from alembic import op
import sqlalchemy as sa


from expression_atlas_db import base, load_db, settings

# revision identifiers, used by Alembic.
revision: str = "c6f2acee9b09"
down_revision: Union[str, None] = "3726b51b1db8"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    base.DataSet.set_alembic(revision)
    if engine_name == "redshift":
        return

    session = base.configure(db_urls["postgres"])()

    studyqueue_table = base.Base.metadata.tables.get("studyqueue")
    studyqueue_table.create(bind=session.bind, checkfirst=True)

    studies = session.query(base.Study).all()
    for s in studies:
        load_db.add_studyqueue(
            s.internal_id,
            session,
            technology="BULK",
            study_id=s.id,
            processed=True,
            status="UPLOADED",
            **{
                c.name: getattr(s, c.name)
                for c in base.StudyQueue.__table__.columns
                if not c.primary_key
                and len(c.foreign_keys) == 0
                and c.name in base.Study.__table__.columns.keys()
                and c.name != "internal_id"
            },
        )
    studyqueues = session.query(base.StudyQueue).all()
    for sq in studyqueues:
        sq.public = True

    session.commit()


def downgrade(engine_name: str, db_urls: Dict[str, str]) -> None:
    if engine_name == "redshift":
        return
    pass
