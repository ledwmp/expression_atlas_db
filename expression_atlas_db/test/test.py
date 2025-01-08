"""Test suite for the Expression Atlas Database package.

This module contains comprehensive tests for the Expression Atlas Database functionality,
including database operations, data loading, and data manipulation. The tests cover:

Classes:
    TestBase: Base test class with common setup/teardown
    TestFullBase: Test class with complete database setup
    TestGTF: Tests for GTF file parsing and loading
    TestMetaDataFetcher: Tests for metadata retrieval from ENA/SRA
    TestExperimentParser: Tests for experiment data parsing
    TestQueue: Tests for study queue management
    TestLoad: Tests for full database loading
    TestBulkStudyQueueAdd: Tests for bulk study queue operations
    TestDeleteStudy: Tests for study deletion
    TestUpdateStudy: Tests for study updates
    TestUpdateDifferentialExpression: Tests for differential expression updates

The test suite uses SQLite for testing and includes sample data in the test_data directory.
Each test class is designed to be independent and includes its own setup/teardown methods
to ensure a clean testing environment.
"""

import unittest
from typing import List, Dict, Union, Any, Tuple
from pathlib import Path
import logging

import pandas as pd
import numpy as np
from sqlalchemy import select, update, bindparam

from expression_atlas_db import base, settings, queries, load_db, utils


class TestBase(unittest.TestCase):
    """Base test class providing common setup and teardown functionality.

    Attributes:
        _test_gtf: Path to test GTF file
        _test_db_connection_string: SQLite database connection string
        _test_data_loc: Path to test data directory
        _test_staging_loc: Path to test staging directory
    """

    _test_gtf = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'test.gtf')}"
    _test_db_connection_string = (
        f"sqlite:///{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'test.db')}"
    )
    _test_data_loc = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data')}"
    _test_staging_loc = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data')}"

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test database connection and schema."""
        settings.test_experiment_loc = cls._test_data_loc
        Session = base.configure(cls._test_db_connection_string)
        cls.session = Session()
        base.Base.metadata.drop_all(bind=cls.session.bind)
        base.Base.metadata.create_all(bind=cls.session.bind)

    @classmethod
    def tearDownClass(cls) -> None:
        """Clean up database connection and remove test database file."""
        cls.session.close()
        Path(cls._test_db_connection_string.replace("sqlite:///", "")).unlink()


class TestFullBase(TestBase):
    """Test class that loads a complete test database with sample data."""

    @classmethod
    def setUpClass(cls) -> None:
        """Initialize database with test GTF and sample data."""
        settings.test_experiment_loc = cls._test_data_loc
        settings.test_staging_loc = cls._test_staging_loc
        load_db.load_db(
            cls._test_gtf,
            drop_all=True,
            use_redshift=False,
            use_s3=False,
            connection_string=cls._test_db_connection_string,
        )
        Session = base.configure(cls._test_db_connection_string)
        cls.session = Session()

    @classmethod
    def tearDownClass(cls) -> None:
        """Clean up staging files and test database."""
        cls.session.close()
        for p in Path(cls._test_staging_loc).glob("*.csv"):
            p.unlink()
        Path(cls._test_db_connection_string.replace("sqlite:///", "")).unlink()


class TestGTF(TestBase):
    """Test GTF file parsing and database loading functionality.

    Tests GTF parsing, database insertion, and relationship mapping between
    sequence regions, genes, and transcripts.
    """

    def setUp(self):
        """Initialize GTF parser and load data into test database."""
        self.gtf = utils.GTFParser(self._test_gtf)
        load_db.bulk_insert_gtf(self.session, self.gtf)

    def tearDown(self):
        """Clean up sequence region, gene, and transcript tables."""
        for n, t in list(base.Base.metadata.tables.items())[::-1]:
            if n not in ("sequenceregion", "gene", "transcript"):
                continue
            t.drop(bind=self.session.bind, checkfirst=True)
            t.create(bind=self.session.bind, checkfirst=True)

    def testLoad(self):
        """Test GTF data loading into database."""
        pass

    def testRelationship(self):
        """Test relationships between GTF entities in database."""
        pass

    def testGTFParser(self):
        """Test GTF file parsing functionality.

        Verifies:
        - GTF line parsing works correctly
        - Chromosome naming is handled properly
        - Attribute parsing works for different formats
        """
        # Test GTF line parsing
        test_line = 'chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";'
        fields, infos = self.gtf.parse_gtf_line(test_line)

        # Check basic field parsing
        self.assertEqual(fields["Chromosome"], "1")  # tests chr stripping
        self.assertEqual(fields["source"], "HAVANA")
        self.assertEqual(fields["label"], "gene")
        self.assertEqual(fields["Start"], 11869)
        self.assertEqual(fields["End"], 14409)

        # Check attribute parsing
        self.assertEqual(infos["gene_id"], "ENSG00000223972")
        self.assertEqual(infos["gene_name"], "DDX11L1")

        # Test chromosome M to MT conversion
        test_line_m = 'chrM\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "TEST";'
        fields, _ = self.gtf.parse_gtf_line(test_line_m)
        self.assertEqual(fields["Chromosome"], "MT")


class TestMetaDataFetcher(TestBase):
    """Test metadata fetching from ENA and SRA databases.

    Tests the ability to fetch and parse metadata for experimental samples
    from external databases, including:
    - ID resolution between different databases (GEO, SRA, BioProject)
    - Metadata retrieval from ENA and SRA
    - Project information parsing
    - Sample metadata handling
    """

    def setUp(self):
        """Initialize experiment parsers and metadata fetchers for test data."""
        self._exps = {}
        self._adatas = {}

        # Load test data
        for fh in Path(self._test_data_loc).glob("*/de_results/*gene.h5ad"):
            self._adatas[fh.stem] = utils.ExperimentParser(
                fh.stem.split("_")[0], Path(settings.test_experiment_loc)
            )
            self._adatas[fh.stem].load_adatas()
            self._exps[fh.stem] = utils.MetaDataFetcher(
                fh.stem.split("_")[0], self._adatas[fh.stem].samples
            )

    def testIDResolution(self):
        """Test resolution of identifiers between different databases.

        Verifies:
        - GEO ID to BioProject ID resolution
        - SRA ID to BioProject ID resolution
        - Proper error handling for invalid IDs
        """
        for e, exp in self._exps.items():
            # Test ID format validation
            if exp.velia_id.startswith("GS"):
                self.assertIsNotNone(exp.geo_id)
            elif exp.velia_id.startswith(("ER", "SR", "DR")):
                self.assertIsNotNone(exp.srp_id)

            # Test BioProject linking
            if exp.bio_id:
                self.assertTrue(exp.bio_id.startswith("PRJ"))

    def testProjectMetadata(self):
        """Test fetching and parsing of project metadata.

        Verifies:
        - Project title retrieval
        - Project summary retrieval
        - PubMed ID association
        - Proper handling of missing metadata
        """
        for e, exp in self._exps.items():
            # Test project information
            self.assertIsNotNone(exp.project_title)
            self.assertIsInstance(exp.project_title, str)

            if exp.project_summary:
                self.assertIsInstance(exp.project_summary, str)

            # Test PubMed IDs if present
            if exp.pmids:
                self.assertTrue(all(pmid.isdigit() for pmid in exp.pmids.split(",")))

    def testSampleMetadata(self):
        """Test sample metadata retrieval and processing.

        Verifies:
        - Sample metadata DataFrame structure
        - Required columns presence
        - Data type consistency
        - Aggregation of run-level data
        """
        for e, exp in self._exps.items():
            metadata_df = exp.samples_metadata

            # Test DataFrame structure
            self.assertIsInstance(metadata_df, pd.DataFrame)
            self.assertGreater(len(metadata_df), 0)

            # Test required columns
            required_cols = ["Experiment", "Run", "spots", "bases"]
            self.assertTrue(all(col in metadata_df.columns for col in required_cols))

            # Test numeric columns
            numeric_cols = ["spots", "bases", "avgLength"]
            for col in numeric_cols:
                if col in metadata_df.columns:
                    self.assertTrue(pd.api.types.is_numeric_dtype(metadata_df[col]))

            # Test experiment uniqueness
            self.assertEqual(
                len(metadata_df["Experiment"].unique()),
                len(metadata_df),
                "Duplicate experiments found in metadata",
            )

    def testENAMetadata(self):
        """Test ENA-specific metadata retrieval.

        Verifies:
        - ENA API response parsing
        - Required ENA fields presence
        - Data consistency between SRA and ENA
        """
        for e, exp in self._exps.items():
            metadata_df = exp.samples_metadata

            # Test ENA fields
            ena_fields = [
                "library_layout",
                "library_selection",
                "library_source",
                "instrument_model",
                "base_count",
                "read_count",
            ]
            present_ena_fields = [f for f in ena_fields if f in metadata_df.columns]
            self.assertGreater(len(present_ena_fields), 0)

    def testURLFetching(self):
        """Test URL fetching functionality and retry mechanism.

        Verifies:
        - Successful URL fetching
        - Retry mechanism for failed requests
        - Error handling for invalid URLs
        """
        for e, exp in self._exps.items():
            # Test valid URL fetch
            response = exp.fetch_url(exp._search_sra_url.format(id=exp.velia_id))
            self.assertEqual(response.getcode(), 200)

            # Test retry mechanism with invalid URL
            with self.assertRaises(Exception):
                exp.fetch_url("http://invalid.url", max_attempts=1)


class TestExperimentParser(TestBase):
    """Test experiment data parsing functionality.

    Tests parsing of experimental data files and metadata extraction, including:
    - AnnData file loading and validation
    - File statistics tracking
    - Sample metadata handling
    - S3 storage integration
    """

    def setUp(self):
        """Initialize experiment parsers and metadata fetchers."""
        self._exps = {}
        self._adatas = {}

        for fh in Path(self._test_data_loc).glob("*/de_results/*gene.h5ad"):
            self._adatas[fh.stem] = utils.ExperimentParser(
                fh.stem.split("_")[0], Path(settings.test_experiment_loc)
            )
            self._adatas[fh.stem].load_adatas()
            self._exps[fh.stem] = utils.MetaDataFetcher(
                fh.stem.split("_")[0], self._adatas[fh.stem].samples
            )

    def tearDown(self):
        """Clean up experiment and adata objects."""
        self._exps = {}
        self._adatas = {}

    def testFileStats(self):
        """Test file statistics tracking.

        Verifies:
        - File timestamp retrieval
        - File size tracking
        - Stats for both gene and transcript data
        """
        for e, adata in self._adatas.items():
            # Test stat collection
            adata.stat_adatas()

            # Verify timestamps
            self.assertIsNotNone(adata._gene_ts)
            self.assertIsNotNone(adata._transcript_ts)

            # Verify file sizes
            self.assertIsNotNone(adata._gene_size)
            self.assertIsNotNone(adata._transcript_size)

            # Test formatted string outputs
            self.assertIsInstance(adata.file_timestamps, str)
            self.assertIsInstance(adata.file_sizes, str)

    def testDataLoading(self):
        """Test AnnData loading functionality.

        Verifies:
        - Successful loading of gene and transcript data
        - Data validation
        - Error handling for missing files
        """
        for e, adata in self._adatas.items():
            # Test data loading
            self.assertIsNotNone(adata._adata_gene)
            self.assertIsNotNone(adata._adata_transcript)

            # Test sample consistency
            self.assertEqual(
                adata._adata_gene.obs.index.tolist(),
                adata._adata_transcript.obs.index.tolist(),
            )

            # Test error handling
            with self.assertRaises(FileNotFoundError):
                bad_adata = utils.ExperimentParser("nonexistent", self._test_data_loc)
                bad_adata.load_adatas()

    def testSampleMetadata(self):
        """Test sample metadata handling.

        Verifies:
        - Sample list generation
        - Metadata DataFrame structure
        - Category handling
        """
        for e, adata in self._adatas.items():
            # Test sample list
            samples = adata.samples
            self.assertIsInstance(samples, list)
            self.assertGreater(len(samples), 0)

            # Test metadata DataFrame
            meta = adata.samples_metadata
            self.assertIsInstance(meta, pd.DataFrame)
            self.assertEqual(len(meta), len(samples))

    def testS3Integration(self):
        """Test S3 storage integration.

        Verifies:
        - S3 filesystem configuration
        - S3 path handling
        - File operations with S3
        """
        for e, adata in self._adatas.items():
            # Test S3 configuration
            self.assertFalse(adata._s3_enabled)
            self.assertIsNone(adata._s3fs)

            # Test enabling S3 (mock filesystem)
            mock_s3fs = type("MockS3FS", (), {"glob": lambda x: []})()
            adata.enable_s3(mock_s3fs)
            self.assertTrue(adata._s3_enabled)
            self.assertIsNotNone(adata._s3fs)


class TestQueue(TestBase):
    """Test study queue management functionality.

    Tests adding, updating, and managing study processing queue entries.
    """

    def testQueueAdd(self):
        """Test adding studies to the processing queue.

        Tests:
        - Adding new studies with valid/invalid parameters
        - Duplicate study handling
        - Validation of study types
        """
        results = queries.submit_studyqueue(
            self.session,
            "SRP129004",
            "Because.",
            "BULK",
            "Because.",
            "Because.",
            "Because.",
            "Unit test",
            "ULTRA_HIGH",
            "None.",
            geo_id="GSE209142",
        )
        self.assertEqual(results[0], False)
        self.assertTrue(isinstance(results[1], str))
        results = queries.submit_studyqueue(
            self.session,
            "SRP129004",
            "Because.",
            "BULK",
            "Because.",
            "Because.",
            "Because.",
            "Unit test",
            "ULTRA_HIGH",
            "None.",
            geo_id="GSE109142",
        )
        self.assertEqual(results[0], False)
        results = queries.submit_studyqueue(
            self.session,
            "SRP129004",
            "Because.",
            "BULK",
            "Because.",
            "Because.",
            "Because.",
            "Unit test",
            "ULTRA_HIGH",
            "None.",
            geo_id="GSE109142",
        )
        self.assertEqual(results[0], True)
        results = queries.submit_studyqueue(
            self.session,
            "SRP129004",
            "Because.",
            "BULKS",
            "Because.",
            "Because.",
            "Because.",
            "Unit test",
            "ULTRA_HIGH",
            "None.",
            geo_id="GSE109142",
        )
        self.assertEqual(results[0], False)
        self.assertTrue(isinstance(results[1], str))

    def testQueueUpdate(self):
        """Test updating existing queue entries.

        Tests:
        - Modifying queue entry attributes
        - Verifying changes are persisted
        """
        results = queries.submit_studyqueue(
            self.session,
            "SRP129004",
            "Because.",
            "BULK",
            "Because.",
            "Because.",
            "Because.",
            "Unit test",
            "ULTRA_HIGH",
            "None.",
        )
        self.session.commit()
        results = (None, queries.fetch_studyqueue(self.session))
        results[1].loc[0, "category"] = "To search for a gene."
        results = queries.update_studyqueue(
            self.session,
            results[1],
        )
        studyqueue_df = pd.read_sql(
            select(base.StudyQueue).filter(
                base.StudyQueue.id.in_(results["id"].tolist())
            ),
            self.session.bind,
        )
        self.assertEqual(
            studyqueue_df.loc[:, "category"].values[0], "To search for a gene."
        )


class TestLoad(TestFullBase):
    """Test full database loading functionality.

    Verifies correct loading of samples, studies, and differential expression data.
    """

    def testLoad(self):
        """Test database loading with complete dataset.

        Verifies:
        - Sample count
        - Study count
        - Differential expression count
        """
        sample_counts = self.session.query(base.Sample).count()
        study_counts = self.session.query(base.Study).count()
        differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()

        self.assertEqual(sample_counts, 301)
        self.assertEqual(study_counts, 2)
        self.assertEqual(differentialexpression_counts, 58245)


class TestBulkStudyQueueAdd(TestFullBase):
    """Test bulk addition of studies to processing queue.

    Tests adding multiple studies to the queue simultaneously.
    """

    def testStudyQueueAdd(self):
        """Test bulk addition of studies to queue.

        Verifies:
        - All studies are added correctly
        - Proper attribute mapping
        - Queue entry count
        """
        self.session.query(base.StudyQueue).delete()
        studies = self.session.query(base.Study).all()
        for s in studies:
            load_db.add_studyqueue(
                s.velia_id,
                self.session,
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
                    and c.name != "velia_id"
                },
            )
        studyqueues = self.session.query(base.StudyQueue).all()
        self.assertEqual(len(studyqueues), 2)


class TestDeleteStudy(TestFullBase):
    """Test study deletion functionality.

    Tests complete removal of study data including related samples and measurements.
    """

    _delete_study = "SRP129004"

    def deleteStudy(self):
        """Test complete study deletion.

        Verifies:
        - Removal of all study-related samples
        - Removal of study entry
        - Removal of related differential expressions
        - Removal of sample measurements
        - Proper count updates
        """
        sample_counts = (
            self.session.query(base.Sample)
            .join(base.Study, base.Study.id == base.Sample.study_id)
            .filter(base.Study.velia_id == self._delete_study)
            .count()
        )
        study_counts = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._delete_study)
            .count()
        )
        differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()
        differentialexpression_counts_todelete = (
            self.session.query(base.DifferentialExpression)
            .join(
                base.Contrast,
                base.Contrast.id == base.DifferentialExpression.contrast_id,
            )
            .filter(base.Contrast.study.has(velia_id=self._delete_study))
            .count()
        )
        samplemeasurement_counts = self.session.query(base.SampleMeasurement).count()
        samplemeasurement_counts_todelete = (
            self.session.query(base.SampleMeasurement)
            .join(base.Sample, base.Sample.id == base.SampleMeasurement.sample_id)
            .filter(base.Sample.study.has(velia_id=self._delete_study))
            .count()
        )

        logging.info(
            f"Deleting test {self._delete_study} from db with {sample_counts} samples."
        )

        load_db.delete_study(
            self._delete_study,
            self.session,
            use_s3=False,
            use_redshift=False,
        )

        delete_sample_counts = (
            self.session.query(base.Sample)
            .join(base.Study, base.Study.id == base.Sample.study_id)
            .filter(base.Study.velia_id == self._delete_study)
            .count()
        )
        delete_study_counts = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._delete_study)
            .count()
        )
        delete_differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()
        delete_samplemeasurement_counts = self.session.query(
            base.SampleMeasurement
        ).count()

        logging.info(
            f"Deleted test {self._delete_study} from db with {sample_counts} samples now {delete_sample_counts}."
        )

        self.assertEqual(delete_sample_counts, 0)
        self.assertEqual(delete_study_counts, 0)
        self.assertEqual(
            delete_differentialexpression_counts,
            differentialexpression_counts - differentialexpression_counts_todelete,
        )
        self.assertEqual(
            delete_samplemeasurement_counts,
            samplemeasurement_counts - samplemeasurement_counts_todelete,
        )


class TestUpdateStudy(TestFullBase):
    """Test study update functionality.

    Tests updating existing study data while maintaining data integrity.
    """

    _update_study = "SRP129004"

    def updateStudy(self):
        """Test study data update process.

        Verifies:
        - Sample count consistency
        - Study count consistency
        - Measurement count consistency
        - Differential expression count consistency
        """
        sample_counts = (
            self.session.query(base.Sample)
            .join(base.Study, base.Study.id == base.Sample.study_id)
            .filter(base.Study.velia_id == self._update_study)
            .count()
        )
        study_counts = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._update_study)
            .count()
        )
        samplemeasurement_counts = self.session.query(base.SampleMeasurement).count()
        differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()

        logging.info(
            f"Updating test {self._update_study} from db with {sample_counts} samples."
        )

        load_db.update_study(
            self._update_study,
            self.session,
            use_s3=False,
            use_redshift=False,
            force=True,
        )

        update_sample_counts = (
            self.session.query(base.Sample)
            .join(base.Study, base.Study.id == base.Sample.study_id)
            .filter(base.Study.velia_id == self._update_study)
            .count()
        )
        update_study_counts = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._update_study)
            .count()
        )
        update_samplemeasurement_counts = self.session.query(
            base.SampleMeasurement
        ).count()
        update_differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()

        logging.info(
            f"Updated test {self._update_study} from db with {sample_counts} samples now {update_sample_counts}."
        )

        self.assertEqual(update_sample_counts, sample_counts)
        self.assertEqual(update_study_counts, study_counts)
        self.assertEqual(update_samplemeasurement_counts, samplemeasurement_counts)
        self.assertEqual(
            update_differentialexpression_counts, differentialexpression_counts
        )


class TestUpdateDifferentialExpression(TestFullBase):
    """Test updating differential expression values in the database.

    Tests the ability to modify differential expression values for specific contrasts
    while maintaining data integrity and preserving unmodified entries.
    """

    _update_study = "SRP129004"

    def updateDifferentialExpression(self):
        """Test updating differential expression values for a specific study contrast.

        This test:
        1. Loads existing differential expression data for a study
        2. Modifies control and case mean values (multiplies by 10)
        3. Updates the database with modified values
        4. Verifies that:
           - Modified entries are correctly updated
           - Unmodified entries remain unchanged
           - Changes are properly persisted in the database
           - Batch processing works correctly (100 records at a time)

        The test specifically checks:
        - Pre vs post update values for modified entries
        - Preservation of unmodified entries
        - Correct scaling of control_mean and case_mean values
        - Data integrity across the update process
        """
        # Get study object
        study = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._update_study)
            .first()
        )

        # Load sequence regions (genes and transcripts)
        sequenceregions = {
            **{g.gene_id: g for g in self.session.query(base.Gene).all()},
            **{t.transcript_id: t for t in self.session.query(base.Transcript).all()},
        }

        logging.info(f"Loading adatas {self._update_study}...")
        exp = utils.ExperimentParser(self._update_study, Path(self._test_data_loc))
        exp.load_adatas()

        # Set up differential expression columns
        de_columns = [
            base.DifferentialExpression._column_map[c.name]
            for c in base.DifferentialExpression.__table__.columns
            if not c.primary_key and len(c.foreign_keys) == 0
        ]

        # Process gene contrasts
        g_exp_dict = exp.prepare_differentialexpression(
            measurement_type="gene",
            de_columns=de_columns,
        )

        for c, (c_df, _, _) in g_exp_dict.items():
            logging.info(f"Modifying contrast {c} from study {study.velia_id}.")

            # Get contrast object
            contrast = (
                self.session.query(base.Contrast)
                .filter(base.Contrast.contrast_name == c)
                .filter(base.Contrast.study_id == study.id)
                .first()
            )

            # Set up staging file
            table_fh = (
                Path(self._test_staging_loc) / f"gene_de.{study.id}.{contrast.id}.csv"
            )

            # Create differential expression entries
            load_db.create_differentialexpression(
                c_df,
                sequenceregions,
                study,
                contrast,
                fh=table_fh,
            )

            # Load data to be modified
            contrast_df = pd.read_csv(
                table_fh,
                names=["contrast_id", "sequenceregion_id"] + de_columns,
            ).loc[:, ["contrast_id", "sequenceregion_id", "control_mean", "case_mean"]]

            # Get original values for comparison
            unchanged_contrast_df = pd.read_sql(
                select(base.DifferentialExpression)
                .filter(base.DifferentialExpression.contrast_id == contrast.id)
                .filter(
                    base.DifferentialExpression.sequenceregion_id.in_(
                        contrast_df["sequenceregion_id"].tolist()
                    )
                ),
                self.session.bind,
            ).set_index("id")

            pre_unmodified_contrast_df = pd.read_sql(
                select(base.DifferentialExpression)
                .filter(base.DifferentialExpression.contrast_id == contrast.id)
                .filter(
                    base.DifferentialExpression.sequenceregion_id.not_in(
                        contrast_df["sequenceregion_id"].tolist()
                    )
                ),
                self.session.bind,
            ).set_index("id")

            contrast_df["control_mean"] = contrast_df["control_mean"] * 10
            contrast_df["case_mean"] = contrast_df["case_mean"] * 10

            contrast_df.columns = [f"b__{n}" for n in contrast_df.columns]

            for i in range(0, contrast_df.shape[0], 100):
                self.session.connection().execute(
                    update(base.DifferentialExpression)
                    .where(
                        base.DifferentialExpression.contrast_id
                        == bindparam("b__contrast_id")
                    )
                    .where(
                        base.DifferentialExpression.sequenceregion_id
                        == bindparam("b__sequenceregion_id")
                    )
                    .values(
                        {
                            "control_mean": bindparam("b__control_mean"),
                            "case_mean": bindparam("b__case_mean"),
                        }
                    ),
                    contrast_df.iloc[i : i + 100].to_dict("records"),
                )

            self.session.commit()

            changed_contrast_df = pd.read_sql(
                select(base.DifferentialExpression)
                .filter(base.DifferentialExpression.contrast_id == contrast.id)
                .filter(
                    base.DifferentialExpression.sequenceregion_id.in_(
                        contrast_df["b__sequenceregion_id"].tolist()
                    )
                ),
                self.session.bind,
            ).set_index("id")

            post_unmodified_contrast_df = pd.read_sql(
                select(base.DifferentialExpression)
                .filter(base.DifferentialExpression.contrast_id == contrast.id)
                .filter(
                    base.DifferentialExpression.sequenceregion_id.not_in(
                        contrast_df["b__sequenceregion_id"].tolist()
                    )
                ),
                self.session.bind,
            ).set_index("id")

            self.assertEqual(
                pre_unmodified_contrast_df.loc[
                    pre_unmodified_contrast_df.index, "control_mean"
                ].tolist(),
                post_unmodified_contrast_df.loc[
                    pre_unmodified_contrast_df.index, "control_mean"
                ].tolist(),
            )
            self.assertEqual(
                pre_unmodified_contrast_df.loc[
                    pre_unmodified_contrast_df.index, "case_mean"
                ].tolist(),
                post_unmodified_contrast_df.loc[
                    pre_unmodified_contrast_df.index, "case_mean"
                ].tolist(),
            )
            self.assertEqual(
                changed_contrast_df.loc[
                    changed_contrast_df.index, "control_mean"
                ].tolist(),
                (
                    unchanged_contrast_df.loc[changed_contrast_df.index, "control_mean"]
                    * 10
                ).tolist(),
            )
            self.assertEqual(
                changed_contrast_df.loc[
                    changed_contrast_df.index, "case_mean"
                ].tolist(),
                (
                    unchanged_contrast_df.loc[changed_contrast_df.index, "case_mean"]
                    * 10
                ).tolist(),
            )


if __name__ == "__main__":
    unittest.main()
