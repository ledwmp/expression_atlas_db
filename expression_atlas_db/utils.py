"""Utility classes for parsing and fetching biological data from various sources.

This module provides classes for handling different types of biological data files and metadata:

Classes:
    GTFParser: Parser for GTF (Gene Transfer Format) files that extracts gene and 
               transcript information into structured DataFrames.
    
    ExperimentParser: Handles loading and processing of gene/transcript expression data 
                     from AnnData files, including differential expression results.
    
    MetaDataFetcher: Retrieves and parses metadata from bioinformatics databases 
                    (NCBI SRA, BioProject, ENA) using their respective APIs.

The module is designed to work with the Expression Atlas database system and supports
both local filesystem and S3 storage access for data files.
"""

import sys
from typing import List, Tuple, Dict, Union, Optional
from pathlib import Path
from io import StringIO
import time
import http
import logging
import tempfile
import urllib
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import s3fs
import anndata as ad

from expression_atlas_db import settings


class GTFParser:
    """Parser for GTF (Gene Transfer Format) files.

    Reads and parses GTF files to extract gene and transcript information into structured
    pandas DataFrames.

    Attributes:
        gtf_path (Path): Path to the GTF file
        gtf_name (str): Name of the GTF file without extension
        genome_build (str | None): Genome build identifier from GTF header
        genome_accession (str | None): Genome accession from GTF header
        transcript_df (pd.DataFrame): DataFrame containing transcript information
        gene_df (pd.DataFrame): DataFrame containing gene information
    """

    def __init__(
        self, gtf_path: Union[str, Path], drop_duplicate_id: bool = False
    ) -> None:
        """Initialize GTFParser with a GTF file path.

        Args:
            gtf_path (Union[str, Path]): Path to the GTF file
            drop_duplicate_id (bool): If True, removes duplicate gene/transcript IDs keeping first occurrence
        """
        self.gtf_path = Path(gtf_path)
        self.gtf_name = self.gtf_path.stem
        self.genome_build = None
        self.genome_accession = None
        self.parse_gtf(self.gtf_path, drop_duplicate_id=drop_duplicate_id)

    def parse_gtf(self, gtf_path: Path, drop_duplicate_id: bool = False) -> None:
        """Read GTF from file and create table structures.

        Args:
            gtf_path (Path): Path to GTF file
            drop_duplicate_id (bool): If True, removes duplicate gene/transcript IDs keeping first occurrence
        """
        transcripts = []
        genes = []

        with open(gtf_path, "r") as f_in:
            transcript_count = 0
            gene_count = 0
            for line in f_in:
                if line.startswith("#"):
                    if line.startswith("#!genome-build "):
                        self.genome_build = line.split(" ")[-1].strip()
                    elif line.startswith("#!genome-build-accession "):
                        self.genome_accession = line.split(" ")[-1].strip()
                    continue

                fields, infos = self.parse_gtf_line(line.strip())

                if fields["label"] in (
                    "gene",
                    "transcript",
                ):
                    fields.update(infos)
                    if fields["label"] == "transcript":
                        transcripts.append(fields)
                        transcript_count += 1

                    elif fields["label"] == "gene":
                        genes.append(fields)
                        gene_count += 1

        if not self.genome_build:
            self.genome_build = gtf_path.stem.split(".gtf")[0]

        transcript_df = pd.DataFrame(
            transcripts,
            columns=[
                "transcript_id",
                "gene_id",
                "veliadb_id",
                "gene_biotype",
            ],
        )
        gene_df = pd.DataFrame(
            genes,
            columns=[
                "gene_id",
                "veliadb_id",
                "gene_biotype",
            ],
        )

        transcript_df["veliadb_id"] = transcript_df["veliadb_id"].astype(int)
        gene_df["veliadb_id"] = gene_df["veliadb_id"].astype(int)
        transcript_df["assembly_id"] = self.genome_build
        gene_df["assembly_id"] = self.genome_build

        if drop_duplicate_id:
            transcript_df.drop_duplicates("transcript_id", inplace=True)
            gene_df.drop_duplicates("gene_id", inplace=True)

        self.transcript_df = transcript_df.copy()
        self.gene_df = gene_df.copy()

    @staticmethod
    def parse_gtf_line(
        gtf_line: str, no_chr: bool = True
    ) -> Tuple[Dict[str, Union[str, int]], Dict[str, str]]:
        """Parse a single line from a GTF file.

        Args:
            gtf_line (str): Raw line from GTF file
            no_chr (bool): If True, removes 'chr' prefix from chromosome names and converts 'M' to 'MT'

        Returns:
            Tuple[Dict[str, Union[str, int]], Dict[str, str]]: Tuple containing:
                - Dictionary of GTF fields (chromosome, source, label, etc.)
                - Dictionary of additional information from the GTF attributes column
        """
        chrom, source, label, start, end, score, strand, frame, info = gtf_line.split(
            "\t"
        )

        if no_chr:
            if chrom.startswith("chr"):
                chrom = chrom[3:]

            if chrom == "M":
                chrom = "MT"

        fields = {
            "Chromosome": chrom,
            "source": source,
            "label": label,
            "Start": int(start),
            "End": int(end),
            "Score": score,
            "Strand": strand,
            "Frame": frame,
        }
        infos = {
            i.strip()
            .split(" ")[0]
            .strip(): i.strip()
            .split(" ")[1]
            .strip()
            .replace('"', "")
            for i in info.strip(";\n").split(";")
        }
        return fields, infos


class ExperimentParser:
    """Parser for experimental data stored in AnnData format.

    Handles loading and processing of gene and transcript expression data from AnnData files,
    including differential expression results and sample measurements.

    Attributes:
        velia_study_id (str): Identifier for the study
        velia_study_loc (Path): Location of AnnData files
        _s3_enabled (bool): Whether S3 access is enabled
        _s3fs (s3fs.core.S3FileSystem | None): S3 filesystem object if enabled
        _transcript_ts (float | None): Timestamp of transcript data
        _gene_ts (float | None): Timestamp of gene data
        _transcript_size (int | None): Size of transcript data file
        _gene_size (int | None): Size of gene data file
        _adata_gene (ad.AnnData | None): Gene expression AnnData object
        _adata_transcript (ad.AnnData | None): Transcript expression AnnData object
    """

    def __init__(
        self, velia_study_id: str, exp_loc: Path = Path(settings.test_experiment_loc)
    ) -> None:
        """Initialize ExperimentParser with study ID and data location.

        Args:
            velia_study_id (str): Study identifier to load
            exp_loc (Path): Base location of AnnData files
        """
        self.velia_study_id = velia_study_id
        self.velia_study_loc = exp_loc / velia_study_id

        self._s3_enabled = False
        self._s3fs = None

        self._transcript_ts = None
        self._gene_ts = None
        self._transcript_size = None
        self._gene_size = None

        self._adata_gene = None
        self._adata_transcript = None

    @property
    def file_timestamps(self) -> str:
        """Get timestamps of gene and transcript files.

        Returns:
            str: Formatted string of gene/transcript timestamps
        """
        return f"{self._gene_ts}/{self._transcript_ts}"

    @property
    def file_sizes(self) -> str:
        """Get sizes of gene and transcript files.

        Returns:
            str: Formatted string of gene/transcript file sizes
        """
        return f"{self._gene_size}/{self._transcript_size}"

    @property
    def samples(self) -> List[str]:
        """Get list of sample identifiers.

        Returns:
            List[str]: List of sample IDs

        Raises:
            AttributeError: If AnnData objects not loaded
            ValueError: If gene/transcript samples don't match
        """
        if not self._adata_transcript or not self._adata_gene:
            raise AttributeError("Set adatas with load_adatas method.")

        if (
            self._adata_gene.obs.index.values != self._adata_transcript.obs.index.values
        ).all():
            raise ValueError(
                "Study gene/transcript adatas have different samples or misordered samples."
            )

        return self._adata_gene.obs.index.to_list()

    @property
    def samples_metadata(self) -> pd.DataFrame:
        """ """
        if not self._adata_transcript or not self._adata_gene:
            raise AttributeError("Set adatas with load_adatas method.")

        if (
            self._adata_gene.obs.index.values != self._adata_transcript.obs.index.values
        ).all():
            raise ValueError(
                "Study gene/transcript adatas have different samples or misordered samples."
            )

        meta = self._adata_gene.obs.copy()
        meta[meta.columns[meta.dtypes == "category"]] = meta[
            meta.columns[meta.dtypes == "category"]
        ].astype(object)

        return meta

    def enable_s3(self, s3fs: s3fs.core.S3FileSystem) -> None:
        """Enable S3 storage access.

        Args:
            s3fs (s3fs.core.S3FileSystem): Configured S3 filesystem object
        """
        self._s3_enabled = True
        self._s3fs = s3fs

    def stat_adatas(self) -> None:
        """Get stats for both gene and transcript AnnData files.

        Calls stat_adata on both gene and transcript data files to update their
        respective timestamp and size attributes.
        """
        self.stat_adata(glob_pattern="*dds_gene*")
        self.stat_adata(glob_pattern="*dds_transcript*")

    def stat_adata(self, glob_pattern: str) -> Union[str, Path]:
        """Get stats for a specific AnnData file.

        Updates the transcript or gene timestamp and size attributes based on the file stats.

        Args:
            glob_pattern (str): Pattern to match either gene or transcript AnnData files

        Returns:
            Union[str, Path]: Path to the matched AnnData file

        Raises:
            FileNotFoundError: If no file matches the glob pattern
        """
        try:
            if not self._s3_enabled:
                fh = list(Path(self.velia_study_loc / "de_results").glob(glob_pattern))[
                    0
                ]
            else:
                s3_path = str(
                    Path(self.velia_study_loc / "de_results" / glob_pattern)
                ).replace("s3:/", "s3://")
                fh = self._s3fs.glob(s3_path)[0]
        except IndexError:
            raise FileNotFoundError(
                f"AnnData with glob pattern: {glob_pattern} not found for study {self.velia_study_loc}."
            )

        if "transcript" in glob_pattern:
            if not self._s3_enabled:
                self._transcript_ts = fh.stat().st_ctime
                self._transcript_size = fh.stat().st_size
            else:
                self._transcript_ts = self._s3fs.info(fh)["LastModified"].timestamp()
                self._transcript_size = self._s3fs.info(fh)["Size"]
        else:
            if not self._s3_enabled:
                self._gene_ts = fh.stat().st_ctime
                self._gene_size = fh.stat().st_size
            else:
                self._gene_ts = self._s3fs.info(fh)["LastModified"].timestamp()
                self._gene_size = self._s3fs.info(fh)["Size"]

        return fh

    def load_adatas(self) -> None:
        """Load both gene and transcript AnnData objects.

        Populates the _adata_gene and _adata_transcript attributes by calling load_adata
        for each data type.
        """
        self._adata_gene = self.load_adata(adata_type="gene")
        self._adata_transcript = self.load_adata(adata_type="transcript")

    def load_adata(self, adata_type: str = "gene") -> ad.AnnData:
        """Load a single AnnData object and update its stats.

        Args:
            adata_type (str): Type of data to load ("gene" or "transcript")

        Returns:
            ad.AnnData: Loaded AnnData object

        Raises:
            FileNotFoundError: If no matching AnnData file is found
            ValueError: If adata_type is invalid
        """
        glob_pattern = f"*dds_{adata_type}*"
        try:
            if not self._s3_enabled:
                fh = list(Path(self.velia_study_loc / "de_results").glob(glob_pattern))[
                    0
                ]
            else:
                s3_path = str(
                    Path(self.velia_study_loc / "de_results" / glob_pattern)
                ).replace("s3:/", "s3://")
                fh = self._s3fs.glob(s3_path)[0]
        except IndexError:
            raise FileNotFoundError(
                f"AnnData with glob pattern: {glob_pattern} not found for study {self.velia_study_loc}."
            )

        if adata_type == "transcript":
            if not self._s3_enabled:
                self._transcript_ts = fh.stat().st_ctime
                self._transcript_size = fh.stat().st_size
            else:
                self._transcript_ts = self._s3fs.info(fh)["LastModified"].timestamp()
                self._transcript_size = self._s3fs.info(fh)["Size"]
        else:
            if not self._s3_enabled:
                self._gene_ts = fh.stat().st_ctime
                self._gene_size = fh.stat().st_size
            else:
                self._gene_ts = self._s3fs.info(fh)["LastModified"].timestamp()
                self._gene_size = self._s3fs.info(fh)["Size"]

        if self._s3_enabled:
            with tempfile.TemporaryDirectory() as tempdir:
                temp_fh = f"{tempdir}/adata.h5ad"
                self._s3fs.get_file(fh, temp_fh)
                adata = ad.read_h5ad(temp_fh)
        else:
            adata = ad.read_h5ad(fh, backed="r")

        return adata

    def prepare_differentialexpression(
        self,
        de_columns: List[str] = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "-log10_padj",
            "control_mean",
            "case_mean",
        ],
        measurement_type: str = "gene",
    ) -> Dict[str, Tuple[pd.DataFrame, List[str], List[str]]]:
        """Extracts and transforms the contrast data in the adatas in preparation for loading into the
        differentialexpression table.

        Args:
            de_columns (List[str]): The columns to pull out of the contrast dataframes in order of appearance in table.
            measurement_type (str): "gene" or "transcript", pull from the respective adata.
        Returns:
            contrast_dict (Dict[str, Tuple[pd.DataFrame, List[str], List[str]]]):
                Dictionary keyed by contrast name and tuple values with the differential expression dataframe, a list of
                samples on the left side of the contrast, and a list of samples on the right side of the contrast.
        """
        if not self._adata_transcript or not self._adata_gene:
            raise AttributeError("Set adatas with load_adatas method.")

        if measurement_type == "transcript":
            adata = self._adata_transcript
        else:
            adata = self._adata_gene

        contrast_dict = {}
        contrasts = [c for c in adata.uns["contrasts"].keys()]
        for c in contrasts:
            use_columns = de_columns.copy()
            if "control_mean" in de_columns and "case_mean" in de_columns:
                try:
                    left_key = f"{adata.uns['contrasts'][c][0]}_{adata.uns['contrasts'][c][1]}_meannormedcounts"
                    right_key = f"{adata.uns['contrasts'][c][0]}_{adata.uns['contrasts'][c][2]}_meannormedcounts"

                    control_loc = np.where(
                        (adata.uns["stat_results"][c].columns == right_key).tolist()
                    )[0][0]
                    case_loc = np.where(
                        (adata.uns["stat_results"][c].columns == left_key).tolist()
                    )[0][0]
                    use_columns[use_columns.index("control_mean")] = adata.uns[
                        "stat_results"
                    ][c].columns[control_loc]
                    use_columns[use_columns.index("case_mean")] = adata.uns[
                        "stat_results"
                    ][c].columns[case_loc]
                except Exception as e:
                    raise Exception("DE columns not formatted properly.") from e
            else:
                use_columns = [
                    c
                    for c in use_columns
                    if c
                    not in (
                        "control_mean",
                        "case_mean",
                    )
                ]
                de_columns = [
                    c
                    for c in use_columns
                    if c
                    not in (
                        "control_mean",
                        "case_mean",
                    )
                ]

            if "log10_pvalue" in de_columns:
                use_columns[use_columns.index("log10_pvalue")] = "pvalue"

            de_df = adata.uns["stat_results"][c][use_columns].copy()

            de_df.columns = [c.lower().strip("-") for c in de_columns]
            de_df["sample_condition_key"] = adata.uns["contrasts"][c][0]
            de_df["left_condition_level"] = adata.uns["contrasts"][c][1]
            de_df["right_condition_level"] = adata.uns["contrasts"][c][2]
            de_df["left_condition_display"] = c.upper().split("_VS_")[0]
            de_df["right_condition_display"] = c.upper().split("_VS_")[1]

            de_df.loc[de_df["padj"].isna(), "padj"] = 1.0

            # TODO: Find a better way to truncate the small pvalues to prevent load errors in redshift.
            # Redshift/Postgres double precision minimum should be ~2.2250738585072014e-308

            de_df.loc[
                np.isneginf(de_df["padj"]) | (de_df["padj"] < sys.float_info.min * 10),
                "padj",
            ] = 0.0
            de_df.loc[
                np.isneginf(de_df["pvalue"])
                | (de_df["pvalue"] < sys.float_info.min * 10),
                "pvalue",
            ] = 0.0
            de_df.loc[de_df["pvalue"].isna(), "pvalue"] = 1.0
            de_df.loc[de_df["log2foldchange"].isna(), "log2foldchange"] = 0.0
            if "-log10_padj" in de_columns:
                de_df["log10_padj"] = -1.0 * de_df["log10_padj"]

            if "log10_pvalue" in de_columns:
                de_df.loc[de_df["log10_pvalue"].isna(), "log10_pvalue"] = 1.0
                de_df["log10_pvalue"] = np.log10(de_df["log10_pvalue"])
                de_df.loc[de_df["log10_pvalue"].isna(), "log10_pvalue"] = 1.0

                # This is the same log10_pvalue limit used in the processing pipeline.

                de_df.loc[de_df["log10_pvalue"] < -400.0, "log10_pvalue"] = -400.0

            left_samples = adata.obs[
                adata.obs[adata.uns["contrasts"][c][0]] == adata.uns["contrasts"][c][1]
            ].index.tolist()
            right_samples = adata.obs[
                adata.obs[adata.uns["contrasts"][c][0]] == adata.uns["contrasts"][c][2]
            ].index.tolist()

            contrast_dict[c.upper()] = (de_df, left_samples, right_samples)

        return contrast_dict

    def prepare_measurements(
        self,
        measurements: List[str] = [
            "counts",
            "normed_counts",
            "raw_tpm",
            "normed_counts_transform",
        ],
        measurement_type: str = "gene",
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """Extract and transform sample measurements from AnnData objects.

        Args:
            measurements (List[str]): List of measurement types to extract from AnnData layers
            measurement_type (str): Type of data to process ("gene" or "transcript")

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]: Tuple containing:
                - Array of sample IDs
                - Array of sequence region IDs
                - 2D array of measurements (samples x measurement types)
                - List of measurement column names

        Raises:
            AttributeError: If AnnData objects not loaded
            ValueError: If measurement_type is invalid
        """
        if not self._adata_transcript or not self._adata_gene:
            raise AttributeError("Set adatas with load_adatas method.")

        if measurement_type == "transcript":
            adata = self._adata_transcript
        else:
            adata = self._adata_gene

        sequenceregions = adata.var.index.values.copy()
        samples = adata.obs.index.values.copy()
        array = np.dstack([adata.layers[l] for l in measurements])

        return (
            np.repeat(samples, sequenceregions.shape[0]),
            np.repeat(
                np.reshape(sequenceregions, (1, -1)), samples.shape[0], axis=0
            ).flatten(),
            np.reshape(array, (-1, len(measurements))),
            measurements,
        )


class MetaDataFetcher:
    """Fetches and parses metadata from various bioinformatics databases.

    Retrieves metadata from NCBI's SRA, BioProject, and ENA databases using their respective APIs.
    Resolves relationships between different database identifiers and extracts study information.

    Attributes:
        _study_id (str): Internal study identifier
        _srx_ids (List[str]): List of SRX experiment identifiers
        _geo_id (str | None): GEO database identifier
        _srp_id (str | None): SRA project identifier
        _project_id (str | None): BioProject identifier
        _project_title (str | None): Project title from BioProject
        _project_summary (str | None): Project summary from BioProject
        _sra_df (pd.DataFrame | None): DataFrame containing SRA metadata
        _pmids (List[str]): List of associated PubMed IDs
    """

    # API URL constants
    _search_sra_url: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term={id}&db=sra"
    )
    _link_bioproject_from_sra_url: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?id={id}&db=bioproject&dbfrom=sra"
    )
    _link_sra_from_bioproject: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?id={id}&db=sra&dbfrom=bioproject"
    )
    _summary_bioproject_url: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?id={id}&db=bioproject"
    )
    _search_bioproject_url: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term={id}&db=bioproject&retmode=xml"
    )
    _fetch_bioproject_url: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id={id}&db=bioproject&retmode=xml"
    )
    _fetch_sra_url_text: str = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id={ids}&db=sra&rettype=runinfo&retmode=text"
    )
    _ena_url: str = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession={bio_id}&fields={extra_params}"
    )

    _ena_fields: List[str] = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "secondary_study_accession",
        "submission_accession",
        "run_alias",
        "experiment_alias",
        "sample_alias",
        "library_layout",
        "library_selection",
        "library_source",
        "library_strategy",
        "library_name",
        "instrument_model",
        "instrument_platform",
        "base_count",
        "read_count",
        "tax_id",
        "scientific_name",
        "sample_title",
        "sample_description",
    ]

    def __init__(self, study_id: str, srx_ids: List[str]) -> None:
        """Initialize MetaDataFetcher with study and experiment identifiers.

        Args:
            study_id (str): Study identifier (GEO or SRA format)
            srx_ids (List[str]): List of SRX experiment identifiers
        """
        self._study_id = study_id
        self._srx_ids = srx_ids
        self._geo_id = None
        self._srp_id = None
        self._project_id = None
        self._project_title = None
        self._project_summary = None
        self._sra_df = None
        self._pmids = []
        self.resolve_all_ids(fetch_srx_info=len(srx_ids) != 0)

    def fetch_url(
        self, url: str, attempt_n: int = 0, max_attempts: int = 5
    ) -> http.client.HTTPResponse:
        """Recursively fetch URLs with retry logic.

        Args:
            url (str): URL to fetch
            attempt_n (int): Current attempt number
            max_attempts (int): Maximum number of retry attempts

        Returns:
            http.client.HTTPResponse: Response from the URL

        Raises:
            Exception: If max attempts reached
            urllib.error.HTTPError: On HTTP error after retries
            urllib.error.URLError: On URL error after retries
        """
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            if "Retry-After" in e.headers:
                time.sleep(int(e.headers["Retry-After"]))
            else:
                time.sleep(15)
            if attempt_n > max_attempts:
                raise Exception(f"Max attempts on url: {url}")
            attempt = attempt_n + 1
            return self.fetch_url(url, attempt_n=attempt)
        except urllib.error.URLError:
            time.sleep(15)
            if attempt_n > max_attempts:
                raise Exception(f"Max attempts on url: {url}")
            attempt = attempt_n + 1
            return self.fetch_url(url, attempt_n=attempt)
        return response

    @property
    def internal_id(self) -> Optional[str]:
        """Get internal study identifier.

        Returns:
            Optional[str]: Study identifier if set
        """
        return self._study_id

    @property
    def geo_id(self) -> Optional[str]:
        """Get GEO database identifier.

        Returns:
            Optional[str]: GEO identifier if available
        """
        return self._geo_id

    @property
    def srp_id(self) -> Optional[str]:
        """Get SRA project identifier.

        Returns:
            Optional[str]: SRA project identifier if available
        """
        return self._srp_id

    @property
    def pmids(self) -> Optional[str]:
        """Get comma-separated list of PubMed IDs.

        Returns:
            Optional[str]: Comma-separated PubMed IDs if available
        """
        return ",".join([pm for pm in self._pmids if pm]) if self._pmids else None

    @property
    def bio_id(self) -> Optional[str]:
        """Get BioProject identifier.

        Returns:
            Optional[str]: BioProject identifier if available
        """
        return self._project_id

    @property
    def project_title(self) -> Optional[str]:
        """Get project title.

        Returns:
            Optional[str]: Project title if available
        """
        return self._project_title

    @property
    def project_summary(self) -> Optional[str]:
        """Get project summary.

        Returns:
            Optional[str]: Project summary if available
        """
        return self._project_summary

    @property
    def samples_metadata(self) -> pd.DataFrame:
        """Get sample metadata DataFrame.

        Returns:
            pd.DataFrame: DataFrame containing sample metadata

        Raises:
            AttributeError: If SRA data not loaded
        """
        if self._sra_df is None:
            raise AttributeError("SRA data not loaded")

        columns_sum = [
            "spots",
            "bases",
            "spots_with_mates",
            "size_MB",
            "base_count",
            "read_count",
        ]
        columns_mean = ["avgLength", "InsertSize", "InsertDev"]

        df = self._sra_df.copy()

        for c in df.columns[~df.columns.isin(columns_sum + columns_mean)]:
            df[c] = df[c].astype(str)

        for c in columns_sum:
            if c not in df.columns:
                continue
            df[c] = df.groupby("Experiment")[c].transform("sum")
        for c in columns_mean:
            if c not in df.columns:
                continue
            df[c] = df.groupby("Experiment")[c].transform("mean")

        df.drop_duplicates("Experiment", inplace=True)

        return df

    def resolve_all_ids(self, fetch_srx_info: bool = True) -> None:
        """Resolve relationships between different database identifiers.

        Links GEO/SRA identifiers to BioProject and fetches associated metadata.

        Args:
            fetch_srx_info (bool): Whether to fetch detailed SRX experiment info

        Raises:
            ValueError: If study_id format is not recognized
        """
        if self._study_id.startswith("GS"):
            self._geo_id = self._study_id
            self._project_id = self.link_project_id(is_sra=False)
        elif any(map(lambda x: self._study_id.startswith(x), ("ER", "SR", "DR"))):
            self._srp_id = self._study_id
            self._project_id = self.link_project_id(is_sra=True)
        else:
            raise ValueError("Cannot define db from study_id.")

        self.fetch_bioproject_info()
        if fetch_srx_info:
            self.fetch_srx_info()

    def link_project_id(self, is_sra: bool = True) -> str:
        """Link study to BioProject identifier.

        Args:
            is_sra (bool): Whether the study ID is from SRA (True) or GEO (False)

        Returns:
            bioproject_id (str): BioProject identifier

        Raises:
            Exception: If unable to link to BioProject
            ValueError: If project ID not found in response
        """
        if is_sra:
            try:
                # Need to get an SRR id from the SRP
                response = self.fetch_url(
                    self._search_sra_url.format(id=self._study_id)
                ).read()
                tree = ET.fromstring(response)
                srr_id = tree.find(".//IdList/Id").text
                response = self.fetch_url(
                    self._link_bioproject_from_sra_url.format(id=srr_id)
                ).read()
                tree = ET.fromstring(response)
                prj_id = tree.find(
                    './/LinkSetDb[LinkName="sra_bioproject_all"]/Link/Id'
                ).text
                if not prj_id:
                    raise ValueError("Unable to find project_id in bioprojects.")
            except Exception as e:
                raise Exception("Unable to link sra_id to bioproject id.") from e
        else:
            try:
                response = self.fetch_url(
                    self._search_bioproject_url.format(id=self._study_id)
                ).read()
                tree = ET.fromstring(response)
                prj_id = tree.find(".//IdList/Id").text
                if not prj_id:
                    raise ValueError("Unable to find project_id in bioprojects.")
            except Exception as e:
                raise Exception("Unable to link geo_id to bioproject id.") from e

        try:
            response = self.fetch_url(
                self._summary_bioproject_url.format(id=prj_id)
            ).read()
            tree = ET.fromstring(response)
            bioproject_id = tree.find(".//Project_Acc").text
            if not bioproject_id:
                raise ValueError(
                    "Unable to find bioproject summary from bioproject id."
                )
        except Exception as e:
            raise Exception(
                "Unable to find bioproject accession from bioproject id."
            ) from e

        return bioproject_id

    def fetch_bioproject_info(self) -> None:
        """Fetch metadata from BioProject database.

        Updates project title, summary, and associated PubMed IDs.

        Raises:
            Exception: If unable to fetch or parse BioProject data
            ValueError: If project ID not found
        """
        try:
            response = self.fetch_url(
                self._search_bioproject_url.format(id=self._project_id)
            ).read()
            tree = ET.fromstring(response)
            project_id = tree.find(".//IdList/Id").text
        except Exception as e:
            raise Exception(
                "Unable to find valid bioproject id given project_id provided."
            ) from e

        try:
            response = self.fetch_url(
                self._fetch_bioproject_url.format(id=project_id)
            ).read()
            tree = ET.fromstring(response)
            if tree.find(".//error"):
                raise ValueError("Cannot find project_id specified.")
        except ValueError as e:
            logging.exception("Unable to fetch bioproject given bioproject id", e)

        try:
            _geo_id = tree.find(".//CenterID").text
            if self._geo_id and self._geo_id != _geo_id:
                raise ValueError(
                    "Provided geo_id not the same as geo_id associated with project."
                )
            elif not self._geo_id:
                self._geo_id = _geo_id
        except Exception as e:
            logging.exception("Unable to link bioproject id to a geo id.", e)

        try:
            self._project_title = tree.find(".//ProjectDescr/Title").text
            self._project_summary = tree.find(".//ProjectDescr/Description").text
        except Exception as e:
            logging.exception(
                "Unable to fetch project name/description from bioproject id.", e
            )

        try:
            for pmid in tree.findall(".//Publication/Reference"):
                self._pmids.append(pmid.text)
        except Exception as e:
            logging.exception(
                "Unable to find pubmed ids attached to bioproject record.", e
            )

    def fetch_srx_info(self, batch_n: int = 50) -> None:
        """Fetch detailed metadata for SRX experiments.

        Retrieves metadata from both SRA and ENA databases and merges the results.

        Args:
            batch_n (int): Number of SRX IDs to process in each batch

        Raises:
            ValueError: If number of retrieved experiments doesn't match input
            Exception: If unable to fetch data from ENA
        """
        # Batch srx_ids into groups of batch_n
        sra_df = pd.DataFrame()
        for i in range(0, len(self._srx_ids), batch_n):
            response = self.fetch_url(
                self._fetch_sra_url_text.format(
                    ids=",".join(self._srx_ids[i : i + batch_n])
                )
            ).read()
            _sra_df = pd.read_csv(
                StringIO(response.decode()), quotechar='"', delimiter=","
            )
            sra_df = pd.concat([sra_df, _sra_df], ignore_index=True)

        if sra_df["Experiment"].unique().shape[0] != len(self._srx_ids):
            raise ValueError(
                "Number of srx_ids retrieved does not equal number initialized. Probably misformatted srx_id."
            )

        if self._study_id.startswith("GS"):
            self._srp_id = sra_df["SRAStudy"].unique()[0]

        try:
            response = self.fetch_url(
                self._ena_url.format(
                    bio_id=self._project_id, extra_params=",".join(self._ena_fields)
                )
            ).read()
            ena_df = pd.read_csv(StringIO(response.decode()), delimiter="\t")
        except Exception as e:
            raise Exception("Unable to fetch data from ena given project id.") from e

        self._sra_df = sra_df.merge(
            ena_df, left_on="Run", right_on="run_accession", how="left"
        )


if __name__ == "__main__":
    pass
