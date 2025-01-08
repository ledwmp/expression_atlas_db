from configparser import ConfigParser
from pathlib import Path
from typing import Dict
from importlib import resources, import_module

config_parser = ConfigParser()

# overwrite defaults settings with settings from the file
with resources.open_text(import_module("expression_atlas_db"), "settings.ini") as f:
    config_parser.read_file(f)
    config = dict(config_parser)

def build_postgres_connection_string(user: str, password: str, host: str, database: str) -> str:
    """Create a PostgreSQL connection string from components."""
    return f"postgresql://{user}:{password}@{host}/{database}"

def build_redshift_connection_string(user: str, password: str, url: str) -> str:
    """Create a Redshift connection string from components."""
    return f"redshift+psycopg2://{user}:{password}@{url}"

# Database connection strings
db_connection_string = build_postgres_connection_string(
    user=config["DATABASE"]["postgres_user"],
    password=config["DATABASE"]["postgres_password"],
    host=config["DATABASE"]["postgres_host"],
    database=config["DATABASE"]["postgres_database"],
)

db_dev_connection_string = build_postgres_connection_string(
    user=config["DATABASE"]["postgres_user"],
    password=config["DATABASE"]["postgres_password"],
    host=config["DATABASE"]["postgres_host"],
    database=config["DATABASE"]["postgres_test_database"],
)

# Redshift connection strings
redshift_connection_string = build_redshift_connection_string(
    user=config["DATABASE"]["redshift_user"],
    password=config["DATABASE"]["redshift_user"], 
    url=config["DATABASE"]["redshift_url"],
)

redshift_dev_connection_string = build_redshift_connection_string(
    user=config["DATABASE"]["redshift_user"],
    password=config["DATABASE"]["redshift_user"],
    url=config["DATABASE"]["redshift_test_url"],
)

# SQLite connection string
sqlite_connection_string = f"sqlite:///{config['DATABASE']['sqlite_test_database']}"

# Data locations
test_experiment_loc = config["DATA"]["test_experiment_loc"]
s3_experiment_loc = config["DATA"]["s3_experiment_loc"]
test_staging_loc = config["DATA"]["test_staging_loc"]
s3_staging_loc = config["DATA"]["s3_staging_loc"]

# AWS configuration
redshift_iam_role = config["DATA"]["redshift_iam_role"]
