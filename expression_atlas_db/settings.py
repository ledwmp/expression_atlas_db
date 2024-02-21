
"""Retrive local user settings"""

from configparser import ConfigParser, NoOptionError
from pathlib import Path
from sys import modules
import os

self = modules[__name__]
settings_ini = (Path(*Path(os.path.abspath(__file__)).parts[:-1]) / 'settings.ini').resolve()

config_parser = ConfigParser()

# overwrite defaults settings with settings from the file
if settings_ini.exists():
    config_parser.read(settings_ini)
    config = dict(config_parser)    
else:
    raise Exception('No settings files at path: %s' % settings_ini)

# set up the database connection string
self.db_connection_string = ('postgresql://%s:%s@%s/%s' %
                             (config['DATABASE']['postgres_user'],
                              config['DATABASE']['postgres_password'],
                              config['DATABASE']['postgres_host'],
                              config['DATABASE']['postgres_database']
                            ))

# set up the redshift connection string 
self.redshift_connection_string = ('redshift+psycopg2://%s:%s@%s' %
                                   (config['DATABASE']['redshift_user'],
                                    config['DATABASE']['redshift_user'],
                                    config['DATABASE']['redshift_url'],
                                   ))

# set test sqlite database connection string 
self.sqlite_connection_string = ('sqlite:///%s' % (config['DATABASE']['sqlite_test_database']))

# set test data location
self.test_experiment_loc = config['DATA']['test_experiment_loc']

# set s3 data location
self.s3_experiment_loc = config['DATA']['s3_experiment_loc']
