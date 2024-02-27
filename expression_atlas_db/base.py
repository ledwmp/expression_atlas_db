from typing import Tuple
from types import MethodType

import sqlalchemy
from sqlalchemy import Column, Integer, String, Float, \
                        Boolean, JSON, Text, ForeignKey, \
                        Identity, Index, MetaData, create_engine
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy.sql.expression import text
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy_redshift.dialect import SUPER

from expression_atlas_db.settings import db_connection_string, \
                                            redshift_connection_string, \
                                            sqlite_connection_string

Base = declarative_base()

class SequenceRegion(Base):
    __tablename__ = 'sequenceregion'

    id = Column(Integer, primary_key=True, autoincrement=True)
    veliadb_id = Column(Integer, nullable=False)
    assembly_id = Column(String(200))
    type = Column(String(30))

    __table_args__ = (UniqueConstraint('veliadb_id','type'),{})
    __mapper_args__ = { 'polymorphic_identity': 'sequence_region', 'polymorphic_on': 'type' }
    
class Gene(SequenceRegion):
    __tablename__ = 'gene'
    
    id = Column(Integer, ForeignKey('sequenceregion.id'), primary_key=True)
    gene_id = Column(String(50))
    gene_biotype = Column(String(200))

    __table_args__ = (UniqueConstraint('gene_id',),{})
    __table_args__ = (Index('idx_gene_id_gene', 'gene_id'),{})
    __mapper_args__ = { 'polymorphic_identity': 'gene' }

class Transcript(SequenceRegion):
    __tablename__ = 'transcript'
    
    id = Column(Integer, ForeignKey('sequenceregion.id'), primary_key=True)
    transcript_id = Column(String(50))
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene_biotype = Column(String(200))
    gene = relationship('Gene', foreign_keys=gene_id)

    __table_args__ = (UniqueConstraint('transcript_id',),{})
    __table_args__ = (Index('idx_transcript_id_transcript', 'transcript_id'),{})
    __mapper_args__ = { 'polymorphic_identity': 'transcript' }

class DataSet(Base):
    __tablename__ = 'dataset'

    id = Column(Integer, primary_key=True, autoincrement=True)
    type = Column(String(30))

    __mapper_args__ = { 'polymorphic_identity': 'dataset', 'polymorphic_on': 'type' }

class Study(DataSet):
    __tablename__ = 'study'
    id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)
    velia_id = Column(String(20), nullable=False)
    geo_id = Column(String(20))
    srp_id = Column(String(20))
    bio_id = Column(String(20))
    pmid = Column(String(100))
    is_public = Column(Boolean, default=True)
    quality = Column(String(20))
    timestamps = Column(String(100))
    sizes = Column(String(100))
    title = Column(String(200))
    description = Column(Text)

    __mapper_args__ = { 'polymorphic_identity': 'study'}

class Sample(DataSet):
    __tablename__ = 'sample'
    
    id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)
    study_id = Column(Integer, ForeignKey('study.id'))
    srx_id = Column(String(20), nullable=False)
    atlas_group = Column(String(200), nullable=False)
    study = relationship('Study', foreign_keys=study_id)
    fields = Column(JSON)

    __table_args__ = (
                    Index('idx_srx_id_sample', 'srx_id'),
                    UniqueConstraint('srx_id','study_id'),
                    {},
                    )
    __mapper_args__ = { 'polymorphic_identity': 'sample' }


class Contrast(DataSet):
    __tablename__ = 'contrast'
    
    id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)
    study_id = Column(Integer, ForeignKey('study.id'))
    contrast_name = Column(String(200), nullable=False)

    # Look into changing these three columns to list types.
    sample_condition_key = Column(String(200))
    left_condition_level = Column(String(200))
    right_condition_level = Column(String(200))

    left_condition_display = Column(String(200))
    right_condition_display = Column(String(200))
    study = relationship('Study', foreign_keys=study_id)


    __table_args__ = (UniqueConstraint('study_id','contrast_name',),{})
    __mapper_args__ = { 'polymorphic_identity': 'contrast' }


class SampleContrast(DataSet):
    __tablename__ = 'samplecontrast'
    
    id = Column(Integer, ForeignKey('dataset.id'), primary_key=True)
    sample_id = Column(Integer, ForeignKey('sample.id'))
    contrast_id = Column(Integer, ForeignKey('contrast.id'))
    study_id = Column(Integer, ForeignKey('study.id'))

    # Look into changing these two levels into list types.
    sample_condition_key = Column(String(200))
    condition_level = Column(String(200))

    condition_display = Column(String(200))
    contrast_side = Column(String(10))
    sample = relationship('Sample', foreign_keys=sample_id, )
    contrast = relationship('Contrast', foreign_keys=contrast_id, )
    study = relationship('Study', foreign_keys=study_id)


    __table_args__ = (UniqueConstraint('sample_id','contrast_id'),{})
    __mapper_args__ = { 'polymorphic_identity': 'samplecontrast' }

class DifferentialExpression(Base):
    __tablename__ = 'differentialexpression'

    id = Column(Integer, primary_key=True, autoincrement=True)
    contrast_id = Column(Integer, ForeignKey('contrast.id'))
    sequenceregion_id = Column(Integer, ForeignKey('sequenceregion.id'))
    basemean = Column(Float)
    log2foldchange = Column(Float)
    lfcse = Column(Float)
    stat = Column(Float)
    pvalue = Column(Float)
    log10_pvalue = Column(Float)
    padj = Column(Float)
    log10_padj = Column(Float)
    control_mean = Column(Float)
    case_mean = Column(Float)
    contrast = relationship('Contrast', foreign_keys=contrast_id, )
    sequenceregion = relationship('SequenceRegion', foreign_keys=sequenceregion_id)

    # __table_args__ = (
    #                 Index('idx_contrast_id_differentialexpression', 'contrast_id'),
    #                 Index('idx_sequenceregion_id_differentialexpression', 'sequenceregion_id'),
    #                 {},
    #                 )
    
    _column_map = {
                'basemean':'baseMean',
                'log2foldchange':'log2FoldChange',
                'lfcse':'lfcSE',
                'stat':'stat',
                'log10_pvalue':'log10_pvalue',
                'pvalue':'pvalue',
                'padj':'padj',
                'log10_padj':'-log10_padj',
                'control_mean':'control_mean',
                'case_mean':'case_mean'
                }

class SampleMeasurement(Base):
    __tablename__ = 'samplemeasurement'

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_id = Column(Integer, ForeignKey('sample.id'))
    sequenceregion_id = Column(Integer, ForeignKey('sequenceregion.id'))
    counts = Column(Float)
    normed_counts = Column(Float)
    tpm = Column(Float)
    normed_counts_transform = Column(Float)
    sample = relationship('Sample', foreign_keys=sample_id)
    sequenceregion = relationship('SequenceRegion', foreign_keys=sequenceregion_id)


    # __table_args__ = (
    #                 Index('idx_sample_id_samplemeasurement', 'sample_id'),
    #                 Index('idx_sequenceregion_id_samplemeasurement', 'sequenceregion_id'),
    #                 {},
    #                 )

    _column_map = {
                'counts':'counts',
                'normed_counts':'normed_counts',
                'tpm':'raw_tpm',
                'normed_counts_transform':'normed_counts_transform'
                }

class _Session(_SA_Session):
    """an sqlalchemy session object to interact with the Velia database
    This object can used to make queries against the database
    The Session will automatically set the search_path to settings.schema
    """

    def __init__(self, *args, **kwargs):
        super(_Session, self).__init__(*args, **kwargs)
        self.get_or_create = MethodType(get_or_create, self)
        self.upsert = MethodType(upsert, self)


    def __repr__(self):
        return "VeliaDB session %d" % (self.__hash__())


def get_or_create(session, class_type, **kwargs):
    """gets an object using filter_by on the unique kwargs. If no such object
    is found in the database, a new one will be created which satisfies
    these constraints. This is why every class that wants to use this
    method to be instantiated needs to have a UniqueConstraint defined.
    """

    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()

    inherited_result = True
    if '__mapper_args__' in class_type.__dict__ and 'inherits' in class_type.__mapper_args__:
        inherited_class_type = class_type.__mapper_args__['inherits']
        for constraint in list(inherited_class_type.__table_args__):
            if constraint.__class__.__name__ == 'UniqueConstraint':
                inherited_unique_cols = constraint.columns.keys()

        try: inherited_result = session.query(inherited_class_type).filter_by(**{k: kwargs[k] for k in inherited_unique_cols}).first()
        except: None

    result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).first()

    if not result or not inherited_result:
        result = class_type(**kwargs)
        session.add(result)
        session.commit()

    return result


def update(session, object, **kwargs):
    """Ideally this would only search on the primary key columns so
    that an update could be made in one call. However, its not currently
    clear how to do that so necessary to pass in the actual object and
    update following a call to get_or_create() There is probably some
    way to do this with class_mapper but its hard right now
    """
    #result = session.query(class_type).filter_by(**kwargs).first()
    #result = session.query(class_type).filter_by(name=kwargs['name']).first()
    #if result is None: return

    for key,value in kwargs.iteritems():
        setattr(object,key,value)
    session.add(object)
    session.commit()

    return object


def upsert(session, class_type, **kwargs):
    """TODO use all unique cols as index
    """
    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()
            
    stmt = insert(class_type).values([kwargs,])
    
    stmt = stmt.on_conflict_do_update(
        index_elements=[eval(f'{class_type.__name__}.{unique_cols[0]}')], set_={k: eval(f'stmt.excluded.{k}', {}, {'stmt': stmt}) for k in kwargs.keys()}
    )
    
    session.execute(stmt)
    session.commit()
    
    result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).first()
    
    return sessionmaker(bind=engine, class_=_Session)

def configure(
            connection_str: str,
            echo:bool=False,
            # ) -> Tuple[sqlalchemy.orm.decl_api.DeclarativeMeta,sqlalchemy.orm.session.Session]:
            ) -> sqlalchemy.orm.session.Session:
    """ 
    """
    engine = create_engine(connection_str, echo=echo)
    # Base.metadata.bind = engine
    Session = sessionmaker(bind=engine, class_=_Session)
    return Session
    
    