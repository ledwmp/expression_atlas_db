from os.path import abspath, dirname
from sys import path
from setuptools import setup, find_packages

setup(
    name='expression_atlas_db',
    version='0.0.1',
    # description="""VeliaDB is a comprehensive resource for internal sORF annotation and experimental data""",
    # author='Stephen Federowicz',
    # author_email='steve@veliatx.com',
    classifiers=[
        'Programming Language :: Python :: 3.8',
    ],
    keywords='microproteins',
    packages=find_packages(),
    install_requires=[
        'SQLAlchemy>=1.3.10,<2.0',
        # 'numpy>=1.17.2',
        'psycopg2-binary',
        # 'scipy>=1.3.1',
        'pytest>=4.6.6',
        'six>=1.12.0',
        'tornado>=4.5.3',
        'configparser>=4.0.2',
        'click',
	'pandas',
    ],
    # entry_points = {
    #     'console_scripts': ['veliadb_load=veliadb.load_db:load_db'],
    # }
)
