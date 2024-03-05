import logging
import re
import time

from sqlalchemy import engine_from_config, pool, create_engine

from alembic import context

from expression_atlas_db import settings, base, load_db

USE_TWOPHASE = False

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

db_names = 'postgres, redshift'

db_urls = {
        'postgres':settings.db_connection_string,
        'redshift':settings.redshift_connection_string,
        }

# target_metadata = {
#             'postgres':base.Base.metadata,
#             'redshift':[t.to_metadata(base.Base.metadata) for n,t in base.Base.metadata.tables.items() \
#                                                 if n in ('samplemeasurement','differentialexpression',)],
#             }
target_metadata = {
            'postgres':'',
            'redshift':'',
            }


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    # for the --sql use case, run migrations for each URL into
    # individual files.

    engines = {}
    for name in re.split(r",\s*", db_names):
        engines[name] = rec = {}
        # rec["url"] = context.config.get_section_option(name, "sqlalchemy.url")
        rec['url'] = db_urls[name]

    for name, rec in engines.items():
        logging.info("Migrating database %s" % name)
        file_ = "%s.sql" % name
        logging.info("Writing output to %s" % file_)
        with open(file_, "w") as buffer:
            context.configure(
                url=rec["url"],
                output_buffer=buffer,
                target_metadata=target_metadata.get(name),
                literal_binds=True,
                dialect_opts={"paramstyle": "named"},
            )
            with context.begin_transaction():
                context.run_migrations(engine_name=name, db_urls=db_urls)


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """

    # for the direct-to-DB use case, start a transaction on all
    # engines, then run all migrations, then commit all transactions.
    engines = {}
    for name in re.split(r",\s*", db_names):
        engines[name] = rec = {}
        # rec["engine"] = engine_from_config(
        #     context.config.get_section(name, {}),
        #     prefix="sqlalchemy.",
        #     poolclass=pool.NullPool,
        # )
        rec["engine"] = create_engine(db_urls[name])

    for name, rec in engines.items():
        engine = rec["engine"]
        rec["connection"] = conn = engine.connect()

        if USE_TWOPHASE:
            rec["transaction"] = conn.begin_twophase()
        else:
            rec["transaction"] = conn.begin()

    try:
        for name, rec in engines.items():
            logging.info("Migrating database %s" % name)
            context.configure(
                connection=rec["connection"],
                upgrade_token="%s_upgrades" % name,
                downgrade_token="%s_downgrades" % name,
                target_metadata=target_metadata.get(name),
            )
            context.run_migrations(engine_name=name, db_urls=db_urls)

        if USE_TWOPHASE:
            for rec in engines.values():
                rec["transaction"].prepare()

        for rec in engines.values():
            rec["transaction"].commit()
    except:
        for rec in engines.values():
            rec["transaction"].rollback()
        raise
    finally:
        for rec in engines.values():
            rec["connection"].close()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
