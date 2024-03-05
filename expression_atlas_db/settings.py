from configparser import ConfigParser
from pathlib import Path
from sys import modules
import os

settings_ini = (
    Path(*Path(os.path.abspath(__file__)).parts[:-2]) / "settings.ini"
).resolve()

config_parser = ConfigParser()

# overwrite defaults settings with settings from the file
if settings_ini.exists():
    config_parser.read(settings_ini)
    config = dict(config_parser)
else:
    raise Exception("No settings files at path: %s" % settings_ini)

# set up the database connection string
db_connection_string = "postgresql://%s:%s@%s/%s" % (
    config["DATABASE"]["postgres_user"],
    config["DATABASE"]["postgres_password"],
    config["DATABASE"]["postgres_host"],
    config["DATABASE"]["postgres_database"],
)

# set up the test database connection string
db_dev_connection_string = "postgresql://%s:%s@%s/%s" % (
    config["DATABASE"]["postgres_user"],
    config["DATABASE"]["postgres_password"],
    config["DATABASE"]["postgres_host"],
    config["DATABASE"]["postgres_test_database"],
)

# set up the redshift connection string
redshift_connection_string = "redshift+psycopg2://%s:%s@%s" % (
    config["DATABASE"]["redshift_user"],
    config["DATABASE"]["redshift_user"],
    config["DATABASE"]["redshift_url"],
)

redshift_dev_connection_string = "redshift+psycopg2://%s:%s@%s" % (
    config["DATABASE"]["redshift_user"],
    config["DATABASE"]["redshift_user"],
    config["DATABASE"]["redshift_test_url"],
)

# set sqlite database connection string
sqlite_connection_string = "sqlite:///%s" % (config["DATABASE"]["sqlite_test_database"])

# set test data location
test_experiment_loc = config["DATA"]["test_experiment_loc"]

# set s3 data location
s3_experiment_loc = config["DATA"]["s3_experiment_loc"]

# set test data staging location
test_staging_loc = config["DATA"]["test_staging_loc"]

# set s3 staging location
s3_staging_loc = config["DATA"]["s3_staging_loc"]

# set redshift iam role
redshift_iam_role = config["DATA"]["redshift_iam_role"]
