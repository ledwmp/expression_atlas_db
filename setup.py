from os.path import abspath, dirname
from sys import path
from setuptools import setup, find_packages

setup(
    name='expression_atlas_db',
    version='0.0.1',
    classifiers=[
        'Programming Language :: Python :: 3.8',
    ],
    packages=find_packages(),
    install_requires=[
        'SQLAlchemy>=1.3.10,<2.0',
        # 'numpy>=1.17.2',
        'psycopg2-binary',
        'pytest>=4.6.6',
        'six>=1.12.0',
        'tornado>=4.5.3',
        'configparser>=4.0.2',
        'click',
	    'pandas',
    ],
)
