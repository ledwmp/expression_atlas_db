import os
import unittest
from typing import List, Dict, Union, Any, Tuple
from pathlib import Path
import logging

import pandas as pd
import numpy as np
import sqlalchemy
from sqlalchemy import select, update, bindparam

from expression_atlas_db import base, settings, queries, load_db, utils


class TestBase(unittest.TestCase):
    """ """

    _test_gtf = (
        f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'veliadb_test.gtf')}"
    )
    _test_db_connection_string = (
        f"sqlite:///{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'test.db')}"
    )
    _test_data_loc = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'test')}"
    _test_staging_loc = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data')}"

    @classmethod
    def setUpClass(cls) -> None:
        """ """
        settings.test_experiment_loc = cls._test_data_loc
        Session = base.configure(cls._test_db_connection_string)
        cls.session = Session()
        base.Base.metadata.drop_all(bind=cls.session.bind)
        base.Base.metadata.create_all(bind=cls.session.bind)

    @classmethod
    def tearDownClass(cls) -> None:
        """ """
        cls.session.close()
        Path(cls._test_db_connection_string.replace("sqlite:///", "")).unlink()


class TestFullBase(TestBase):
    """ """

    @classmethod
    def setUpClass(cls) -> None:
        """ """
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
        """ """
        cls.session.close()
        for p in Path(cls._test_staging_loc).glob("*.csv"):
            p.unlink()
        Path(cls._test_db_connection_string.replace("sqlite:///", "")).unlink()


class TestGTF(TestBase):
    """ """

    _test_v2_gtf = f"{Path('/'.join(Path(__file__).parts[:-1]), 'test_data', 'veliadb_test_v2.gtf')}"

    def setUp(self):
        """ """
        self.gtf = utils.GTFParser(self._test_gtf)
        load_db.bulk_insert_gtf(self.session, self.gtf)

    def tearDown(self):
        """ """
        for n, t in list(base.Base.metadata.tables.items())[::-1]:
            if n not in ("sequenceregion", "gene", "transcript"):
                continue
            t.drop(bind=self.session.bind, checkfirst=True)
            t.create(bind=self.session.bind, checkfirst=True)

    def testLoad(self):
        """ """
        self.assertEqual(self.session.query(base.Gene).count(), 2119)
        self.assertEqual(self.session.query(base.Transcript).count(), 9111)
        self.assertEqual(self.session.query(base.SequenceRegion).count(), 11230)

    def testRelationship(self):
        """ """
        pass

    def testGTFParser(self):
        """ """
        pass

    def testLoadVersionedGTF(self):
        """ """
        v2_gtf = utils.GTFParser(self._test_v2_gtf)
        load_db.bulk_insert_gtf(self.session, v2_gtf)

        self.assertEqual(
            self.session.query(base.Gene)
            .filter(base.Gene.assembly_id == "veliadb_test")
            .count(),
            2119,
        )
        self.assertEqual(
            self.session.query(base.Transcript)
            .filter(base.Transcript.assembly_id == "veliadb_test")
            .count(),
            9111,
        )
        self.assertEqual(
            self.session.query(base.Gene)
            .filter(base.Gene.assembly_id == "veliadb_test_v2")
            .count(),
            2119,
        )
        self.assertEqual(
            self.session.query(base.Transcript)
            .filter(base.Transcript.assembly_id == "veliadb_test_v2")
            .count(),
            9111,
        )

        # Try to re-load the same base gtf, should throw an Integrity Error based on a unique-constraint exception.
        try:
            v1_gtf = utils.GTFParser(self._test_gtf)
            load_db.bulk_insert_gtf(self.session, v1_gtf)
        except sqlalchemy.exc.IntegrityError as e:
            self.session.rollback()
        except Exception as e:
            raise e

        # Test versioned queries against the transcript table.
        all_transcript_srs = queries.fetch_sequenceregions(
            self.session,
            sequenceregions_type="transcript",
        )
        v1_transcript_srs = queries.fetch_sequenceregions(
            self.session,
            sequenceregions_type="transcript",
            assembly_id="veliadb_test",
        )
        v2_transcript_srs = queries.fetch_sequenceregions(
            self.session,
            sequenceregions_type="transcript",
            assembly_id="veliadb_test_v2",
        )
        self.assertEqual(all_transcript_srs.shape[0], 18222)
        self.assertEqual(all_transcript_srs["transcript_id"].nunique(), 9111)
        self.assertEqual(v1_transcript_srs.shape[0], 9111)
        self.assertEqual(v2_transcript_srs.shape[0], 9111)


class TestMetaDataFetcher(TestBase):
    """ """

    def setUp(self):
        """ """
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

    def testENAMeta(self):
        for e, exp in self._exps.items():
            pass

    def testSRAMeta(self):
        for e, exp in self._exps.items():
            pass


class TestExperimentParser(TestBase):

    def setUp(self):
        """ """
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
        """ """
        self._exps = {}
        self._adatas = {}


class TestQueue(TestBase):

    def testQueueAdd(self):
        """ """
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
        """ """
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

    def testLoad(self):
        """ """
        sample_counts = self.session.query(base.Sample).count()
        study_counts = self.session.query(base.Study).count()
        differentialexpression_counts = self.session.query(
            base.DifferentialExpression
        ).count()

        self.assertEqual(sample_counts, 301)
        self.assertEqual(study_counts, 2)
        self.assertEqual(differentialexpression_counts, 58245)


class TestBulkStudyQueueAdd(TestFullBase):

    def testStudyQueueAdd(self):
        """ """
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

    _delete_study = "SRP129004"

    def deleteStudy(self):
        """TODO: There's some flawed logic below. Check on querying sample counts by filtering on study id."""

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

    _update_study = "SRP129004"

    def updateStudy(self):
        """ """

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

    _update_study = "SRP129004"

    def updateDifferentialExpression(self):
        """ """

        study = (
            self.session.query(base.Study)
            .filter(base.Study.velia_id == self._update_study)
            .first()
        )

        sequenceregions = {
            **{g.gene_id: g for g in self.session.query(base.Gene).all()},
            **{t.transcript_id: t for t in self.session.query(base.Transcript).all()},
        }

        logging.info(f"Loading adatas {self._update_study}...")
        exp = utils.ExperimentParser(self._update_study, Path(self._test_data_loc))
        exp.load_adatas()

        # Create column names for differentialexpression table.

        de_columns = [
            base.DifferentialExpression._column_map[c.name]
            for c in base.DifferentialExpression.__table__.columns
            if not c.primary_key and len(c.foreign_keys) == 0
        ]

        # Load gene contrasts.

        g_exp_dict = exp.prepare_differentialexpression(
            measurement_type="gene",
            de_columns=de_columns,
        )

        for c, (c_df, _, _) in g_exp_dict.items():

            logging.info(f"Modifying contrast {c} from study {study.velia_id}.")

            contrast = (
                self.session.query(base.Contrast)
                .filter(base.Contrast.contrast_name == c)
                .filter(base.Contrast.study_id == study.id)
                .first()
            )

            table_fh = (
                Path(self._test_staging_loc) / f"gene_de.{study.id}.{contrast.id}.csv"
            )

            load_db.create_differentialexpression(
                c_df,
                sequenceregions,
                study,
                contrast,
                fh=table_fh,
            )

            # Load the fixed entries to df.

            contrast_df = pd.read_csv(
                table_fh,
                names=["contrast_id", "sequenceregion_id"] + de_columns,
            ).loc[:, ["contrast_id", "sequenceregion_id", "control_mean", "case_mean"]]

            # Load out of differentialexpression table directly.

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
