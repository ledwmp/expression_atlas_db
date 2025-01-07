"""Database models and session configuration for the Expression Atlas database.

This module defines the SQLAlchemy ORM models for the Expression Atlas database,
including sequence regions (genes and transcripts), studies, samples, and expression data.
"""

from sqlalchemy import (
    Column,
    Integer,
    String,
    Float,
    Boolean,
    JSON,
    Text,
    DateTime,
    ForeignKey,
    Index,
    create_engine,
    func,
)
from sqlalchemy.orm import relationship, sessionmaker, session
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base

# String length constraints
MAX_ID_LENGTH = 50
MAX_BIOTYPE_LENGTH = 200
MAX_CONDITION_LENGTH = 200
MAX_ASSEMBLY_ID_LENGTH = 200
MAX_VELIA_ID_LENGTH = 20
MAX_GEO_ID_LENGTH = 20
MAX_SRP_ID_LENGTH = 20
MAX_SRX_ID_LENGTH = 20
MAX_TIMESTAMPS_LENGTH = 100

Base = declarative_base()


class SequenceRegion(Base):
    """Base class for genomic sequence regions.

    Represents a region in the genome with a unique identifier and assembly information.
    This is the parent class for Gene and Transcript.

    Attributes:
        id (int): Primary key
        veliadb_id (int): Reference ID in the Velia database
        assembly_id (str): Genome assembly identifier
        type (str): Type of sequence region (used for polymorphic identity)
    """

    __tablename__ = "sequenceregion"

    id = Column(Integer, primary_key=True, autoincrement=True)
    veliadb_id = Column(Integer, nullable=False)
    assembly_id = Column(String(MAX_ASSEMBLY_ID_LENGTH))
    type = Column(String(30))

    __table_args__ = (UniqueConstraint("veliadb_id", "type"), {})
    __mapper_args__ = {
        "polymorphic_identity": "sequence_region",
        "polymorphic_on": "type",
    }


class Gene(SequenceRegion):
    """Represents a gene in the genome.

    Extends SequenceRegion to include gene-specific information and maintains
    relationships with associated transcripts.

    Attributes:
        id (int): Primary key referencing SequenceRegion
        gene_id (str): Unique identifier for the gene
        gene_biotype (str): Biological type of the gene
        transcripts (List[Transcript]): Related transcript records
    """

    __tablename__ = "gene"

    id = Column(
        Integer, ForeignKey("sequenceregion.id", ondelete="cascade"), primary_key=True
    )
    gene_id = Column(String(MAX_ID_LENGTH))
    gene_biotype = Column(String(MAX_BIOTYPE_LENGTH))
    transcripts = relationship(
        "Transcript",
        foreign_keys="Transcript.gene_id",
        back_populates="gene",
        cascade="all, delete",
    )

    __table_args__ = (Index("idx_gene_id_gene", "gene_id"), {})
    __mapper_args__ = {"polymorphic_identity": "gene"}


class Transcript(SequenceRegion):
    """Represents a transcript in the genome.

    Extends SequenceRegion to include transcript-specific information and maintains
    a relationship with its parent gene.

    Attributes:
        id (int): Primary key referencing SequenceRegion
        transcript_id (str): Unique identifier for the transcript
        gene_id (int): Foreign key to the parent gene
        gene_biotype (str): Biological type of the gene
        gene (Gene): Relationship to the parent gene record
    """

    __tablename__ = "transcript"

    id = Column(
        Integer, ForeignKey("sequenceregion.id", ondelete="cascade"), primary_key=True
    )
    transcript_id = Column(String(MAX_ID_LENGTH))
    gene_id = Column(Integer, ForeignKey("gene.id"))
    gene_biotype = Column(String(MAX_BIOTYPE_LENGTH))
    gene = relationship("Gene", foreign_keys=gene_id)

    __table_args__ = (Index("idx_transcript_id_transcript", "transcript_id"), {})
    __mapper_args__ = {"polymorphic_identity": "transcript"}


class DataSet(Base):
    """Base class for all dataset-related tables.

    Provides common attributes and functionality for managing different types of datasets
    with versioning support through alembic.

    Attributes:
        id (int): Primary key
        alembic_id (str): Database migration version identifier
        created_at (datetime): Timestamp of record creation
        type (str): Type of dataset (used for polymorphic identity)
    """

    __tablename__ = "dataset"

    alembic_revision_id = None

    @classmethod
    def set_alembic(cls, alembic_revision_id: str):
        """Set the alembic revision ID for version tracking.

        Args:
            alembic_revision_id (str): The revision ID from alembic
        """
        cls.alembic_revision_id = alembic_revision_id

    id = Column(Integer, primary_key=True, autoincrement=True)
    alembic_id = Column(String(30), default=lambda: __class__.alembic_revision_id)
    created_at = Column(DateTime, server_default=func.now())
    type = Column(String(30))

    __mapper_args__ = {"polymorphic_identity": "dataset", "polymorphic_on": "type"}


class StudyQueue(DataSet):
    """Represents a study in the processing queue.

    Contains metadata about studies that are queued for processing, including
    their source identifiers, quality metrics, and processing status.

    Attributes:
        id (int): Primary key referencing DataSet
        velia_id (str): Reference ID in the Velia database
        geo_id (str): GEO database identifier
        srp_id (str): SRA project identifier
        study_id (int): Study identifier
        pmid (str): PubMed ID
        status (str): Current processing status
        quality (str): Quality assessment
        technology (str): Sequencing technology
        title (str): Study title
        description (str): Study description
        category (str): Study category
        requestor (str): Study requestor
        priority (str): Processing priority
        contrast (str): Contrast information
        disease (str): Disease information
        tissue (str): Tissue information
        comments (str): Additional comments
        public (bool): Whether the study is publicly accessible
        processed (bool): Whether processing is complete
    """

    __tablename__ = "studyqueue"
    id = Column(Integer, ForeignKey("dataset.id", ondelete="cascade"), primary_key=True)
    velia_id = Column(String(MAX_VELIA_ID_LENGTH), nullable=False)
    geo_id = Column(String(MAX_GEO_ID_LENGTH))
    srp_id = Column(String(MAX_SRP_ID_LENGTH))
    study_id = Column(Integer)
    pmid = Column(Text)
    status = Column(Text)
    quality = Column(Text)
    technology = Column(Text)
    title = Column(Text)
    description = Column(Text)
    category = Column(Text)
    requestor = Column(Text)
    priority = Column(Text)
    contrast = Column(Text)
    disease = Column(Text)
    tissue = Column(Text)
    comments = Column(Text)
    public = Column(Boolean, default=False)
    processed = Column(Boolean, default=False)

    __mapper_args__ = {"polymorphic_identity": "studyqueue"}


class Study(DataSet):
    """Represents a processed study in the database.

    Contains comprehensive study information including metadata, publication info,
    and relationships to samples, contrasts, and sample contrasts.

    Attributes:
        id (int): Primary key referencing DataSet
        velia_id (str): Reference ID in the Velia database
        geo_id (str): GEO database identifier
        srp_id (str): SRA project identifier
        bio_id (str): Biological identifier
        pmid (str): PubMed ID
        public (bool): Whether study is publicly accessible
        quality (str): Quality assessment
        timestamps (str): Processing timestamps
        sizes (str): Data sizes
        title (str): Study title
        description (str): Study description
        samples (List[Sample]): Related sample records
        contrasts (List[Contrast]): Related contrast records
        samplecontrasts (List[SampleContrast]): Related sample-contrast mapping records
    """

    __tablename__ = "study"
    id = Column(Integer, ForeignKey("dataset.id", ondelete="cascade"), primary_key=True)
    velia_id = Column(String(MAX_VELIA_ID_LENGTH), nullable=False)
    geo_id = Column(String(MAX_GEO_ID_LENGTH))
    srp_id = Column(String(MAX_SRP_ID_LENGTH))
    bio_id = Column(String(MAX_ID_LENGTH))
    pmid = Column(Text)
    public = Column(Boolean, default=False)
    quality = Column(Text)
    timestamps = Column(String(MAX_TIMESTAMPS_LENGTH))
    sizes = Column(String(100))
    title = Column(Text)
    description = Column(Text)
    samples = relationship(
        "Sample",
        foreign_keys="Sample.study_id",
        back_populates="study",
        cascade="all, delete",
    )
    contrasts = relationship(
        "Contrast",
        foreign_keys="Contrast.study_id",
        back_populates="study",
        cascade="all, delete",
    )
    samplecontrasts = relationship(
        "SampleContrast",
        foreign_keys="SampleContrast.study_id",
        back_populates="study",
        cascade="all, delete",
    )

    __mapper_args__ = {"polymorphic_identity": "study"}


class Sample(DataSet):
    """Represents a biological sample in a study.

    Contains sample-specific information including identifiers, grouping information,
    and custom fields for additional metadata.

    Attributes:
        id (int): Primary key referencing DataSet
        study_id (int): Foreign key to the parent study
        srx_id (str): SRA experiment identifier
        atlas_group (str): Sample grouping identifier
        fields (dict): JSON field containing additional metadata
        study (Study): Parent study relationship
        samplecontrasts (List[SampleContrast]): Related sample-contrast mappings
    """

    __tablename__ = "sample"

    id = Column(Integer, ForeignKey("dataset.id", ondelete="cascade"), primary_key=True)
    study_id = Column(Integer, ForeignKey("study.id"))
    srx_id = Column(String(MAX_SRX_ID_LENGTH), nullable=False)
    atlas_group = Column(String(MAX_CONDITION_LENGTH), nullable=False)
    fields = Column(JSON)
    study = relationship("Study", foreign_keys=study_id)
    samplecontrasts = relationship(
        "SampleContrast",
        foreign_keys="SampleContrast.sample_id",
        cascade="all, delete",
        back_populates="sample",
    )

    __table_args__ = (
        Index("idx_srx_id_sample", "srx_id"),
        UniqueConstraint("srx_id", "study_id"),
        {},
    )
    __mapper_args__ = {"polymorphic_identity": "sample"}


class Contrast(DataSet):
    """Represents a comparison between two conditions in a study.

    Defines the parameters for differential expression analysis between
    two specific conditions within a study.

    Attributes:
        id (int): Primary key referencing DataSet
        study_id (int): Foreign key to the parent study
        contrast_name (str): Name identifying the contrast
        sample_condition_key (str): Key used to group samples
        left_condition_level (str): Control condition identifier
        right_condition_level (str): Test condition identifier
        left_condition_display (str): Display name for control condition
        right_condition_display (str): Display name for test condition
        study (Study): Parent study relationship
        samplecontrasts (List[SampleContrast]): Related sample-contrast mappings
    """

    __tablename__ = "contrast"

    id = Column(Integer, ForeignKey("dataset.id", ondelete="cascade"), primary_key=True)
    study_id = Column(Integer, ForeignKey("study.id"))
    contrast_name = Column(String(MAX_CONDITION_LENGTH), nullable=False)
    sample_condition_key = Column(String(MAX_CONDITION_LENGTH))
    left_condition_level = Column(String(MAX_CONDITION_LENGTH))
    right_condition_level = Column(String(MAX_CONDITION_LENGTH))
    left_condition_display = Column(String(MAX_CONDITION_LENGTH))
    right_condition_display = Column(String(MAX_CONDITION_LENGTH))
    study = relationship("Study", foreign_keys=study_id)
    samplecontrasts = relationship(
        "SampleContrast",
        foreign_keys="SampleContrast.contrast_id",
        cascade="all, delete",
        back_populates="contrast",
    )

    __table_args__ = (
        UniqueConstraint(
            "study_id",
            "contrast_name",
        ),
        {},
    )
    __mapper_args__ = {"polymorphic_identity": "contrast"}


class SampleContrast(DataSet):
    """Represents the association between samples and contrasts.

    Maps samples to their respective contrasts and specifies which condition
    (case/control) the sample represents in the contrast.

    Attributes:
        id (int): Primary key referencing DataSet
        sample_id (int): Foreign key to the sample
        contrast_id (int): Foreign key to the contrast
        study_id (int): Foreign key to the study
        sample_condition_key (str): Sample condition grouping key
        condition_level (str): Sample condition identifier
        condition_display (str): Display name for condition
        contrast_side (str): Whether sample is case or control
        sample (Sample): Related sample record
        contrast (Contrast): Related contrast record
        study (Study): Related study record
    """

    __tablename__ = "samplecontrast"

    id = Column(Integer, ForeignKey("dataset.id", ondelete="cascade"), primary_key=True)
    sample_id = Column(Integer, ForeignKey("sample.id"))
    contrast_id = Column(Integer, ForeignKey("contrast.id"))
    study_id = Column(Integer, ForeignKey("study.id"))
    sample_condition_key = Column(String(MAX_CONDITION_LENGTH))
    condition_level = Column(String(MAX_CONDITION_LENGTH))
    condition_display = Column(String(MAX_CONDITION_LENGTH))
    contrast_side = Column(String(10))
    sample = relationship(
        "Sample",
        foreign_keys=sample_id,
    )
    contrast = relationship(
        "Contrast",
        foreign_keys=contrast_id,
    )
    study = relationship("Study", foreign_keys=study_id)

    __table_args__ = (UniqueConstraint("sample_id", "contrast_id"), {})
    __mapper_args__ = {"polymorphic_identity": "samplecontrast"}


class DifferentialExpression(Base):
    """Stores differential expression analysis results.

    Contains statistical measures and expression changes for genes/transcripts
    between conditions in a contrast.

    Attributes:
        id (int): Primary key
        contrast_id (int): Foreign key to the contrast
        sequenceregion_id (int): Foreign key to the gene/transcript
        basemean (float): Mean expression level
        log2foldchange (float): Log2 fold change between conditions
        lfcse (float): Log fold change standard error
        stat (float): Test statistic
        pvalue (float): Statistical significance
        log10_pvalue (float): -Log10 transformed p-value
        padj (float): Adjusted p-value
        log10_padj (float): -Log10 transformed adjusted p-value
        control_mean (float): Mean expression in control condition
        case_mean (float): Mean expression in case condition
        contrast (Contrast): Related contrast record
        sequenceregion (SequenceRegion): Related sequence region record
    """

    __tablename__ = "differentialexpression"

    id = Column(Integer, primary_key=True, autoincrement=True)
    contrast_id = Column(Integer, ForeignKey("contrast.id"))
    sequenceregion_id = Column(Integer, ForeignKey("sequenceregion.id"))
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
    contrast = relationship(
        "Contrast",
        foreign_keys=contrast_id,
    )
    sequenceregion = relationship("SequenceRegion", foreign_keys=sequenceregion_id)

    _column_map = {
        "basemean": "baseMean",
        "log2foldchange": "log2FoldChange",
        "lfcse": "lfcSE",
        "stat": "stat",
        "log10_pvalue": "log10_pvalue",
        "pvalue": "pvalue",
        "padj": "padj",
        "log10_padj": "-log10_padj",
        "control_mean": "control_mean",
        "case_mean": "case_mean",
    }


class SampleMeasurement(Base):
    """Stores expression measurements for individual samples.

    Contains various normalized and raw expression measurements for
    genes/transcripts in specific samples.

    Attributes:
        id (int): Primary key
        sample_id (int): Foreign key to the sample
        sequenceregion_id (int): Foreign key to the gene/transcript
        counts (float): Raw read counts
        normed_counts (float): Normalized read counts
        tpm (float): Transcripts Per Million
        normed_counts_transform (float): Transformed normalized counts
        sample (Sample): Related sample record
        sequenceregion (SequenceRegion): Related sequence region record
    """

    __tablename__ = "samplemeasurement"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_id = Column(Integer, ForeignKey("sample.id"))
    sequenceregion_id = Column(Integer, ForeignKey("sequenceregion.id"))
    counts = Column(Float)
    normed_counts = Column(Float)
    tpm = Column(Float)
    normed_counts_transform = Column(Float)
    sample = relationship("Sample", foreign_keys=sample_id)
    sequenceregion = relationship("SequenceRegion", foreign_keys=sequenceregion_id)

    _column_map = {
        "counts": "counts",
        "normed_counts": "normed_counts",
        "tpm": "raw_tpm",
        "normed_counts_transform": "normed_counts_transform",
    }


class _Session(session.Session):
    """Custom session class for Expression Atlas database connections.

    Extends SQLAlchemy's Session class with Expression Atlas specific functionality.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return "ExpressionAtlasDB session %d" % (self.__hash__())


def configure(
    connection_str: str,
    echo: bool = False,
    read_only: bool = False,
) -> _Session:
    """Configure and create a database session.

    Args:
        connection_str (str): Database connection string
        echo (bool): Enable SQL query logging. Defaults to False.
        read_only (bool): Set connection to read-only mode. Defaults to False.

    Returns:
        Session (_Session): Configured database session
    """
    engine = create_engine(
        connection_str,
        echo=echo,
        execution_options={"postgresql_readonly": read_only},
        pool_pre_ping=True,
    )
    Session = sessionmaker(bind=engine, class_=_Session)
    return Session
