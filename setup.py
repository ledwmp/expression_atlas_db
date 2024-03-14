from setuptools import setup, find_packages

setup(
    name="expression_atlas_db",
    version="0.0.1",
    classifiers=[
        "Programming Language :: Python :: 3.8",
    ],
    packages=find_packages(),
    install_requires=[
        "SQLAlchemy>=1.3.10,<2.0",
        'numpy>=1.24.4',
        "psycopg2-binary",
        "six>=1.12.0",
        "configparser>=4.0.2",
        "pandas",
        "alembic",
        "sqlalchemy_redshift",
        "s3fs",
        "anndata==0.9.2",
        "defusedxml",
    ],
)
