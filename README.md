# Expression Atlas Database

A database system for managing and analyzing RNA sequencing expression data with integrated differential expression results.

## Overview
The Expression Atlas Database provides a centralized system for storing and accessing RNA-seq expression data and differential expression analysis results. It integrates with veliadb for gene/transcript annotations and automatically fetches metadata from ENA (European Nucleotide Archive) and SRA (Sequence Read Archive).

## Requirements
- AWS S3 access
- AWS Redshift access for large tables
- Postgres or sqlite databse for metadata tables
- Python environment with alembic
- Access to veliadb

## Database Structure
The database consists of several key tables:

### Core Tables
- `studies` - Core study information, QC status, and metadata
- `studyqueue` - Pipeline processing queue and status tracking
- `dataset` - Base table for all dataset-related tables with version tracking

### Sequence Tables
- `sequenceregion` - Base table for genomic regions
- `gene` - Gene records with biotype information
- `transcript` - Transcript records linked to genes

### Expression Data Tables
- `samplemeasurement` - Expression measurements (counts, TPM, etc.) for each sample-feature pair
- `differentialexpression` - Statistical results from differential expression analysis

### Study Organization Tables
- `sample` - Sample information and metadata
- `contrast` - Comparison definitions between conditions
- `samplecontrast` - Maps samples to contrasts and defines case/control relationships

Each table inherits from either `Base` or `DataSet` base classes, providing common functionality like version tracking and timestamps. The database uses foreign key relationships to maintain data integrity and enable efficient querying across the expression atlas.

## Usage

### Study Ingestion
Studies can be ingested into the database through the following process:

1. **Data Preparation**
   - Prepare study data using the `de_processing` notebooks (01-04) in the expression_atlas repo
   - Ensure data meets the required format specifications

2. **S3 Upload**
   Copy study directory to S3 at: `<s3_experiment_loc>` with structure:
   ```
   SRP000001  
   |   ...
   └───fetchngs_output
   │   │   ...
   └───rnaseq_output
   │   │   ...
   └───de_results
       │   SRP000001_dds_transcript.h5_ad
       │   SRP000001_dds_gene.h5_ad
       |   ...
   ``` 

3. **Data Requirements**
   AnnData files must meet these specifications:
   - obs indexed on srx accession id
   - var indexed on gene_id/transcript_id from veliadb gtf
   - Required layers: `["counts", "normed_counts", "raw_tpm", "normed_counts_transform"]`
   - uns dictionary requirements:
     - "contrasts" with format: `[<factor name>, <test factor>, <reference factor>]`
     - "stat_results" with specified column structure: `["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
        "padj", "-log10_padj"]` and an additional two columns with case and 
        control means of normed_counts_transform named like: `<factor 
        name>_<factor level>_meannormedcounts`

4. **Database Loading**
   Run from project root:
   ```bash
   alembic revision -m 'add'
   alembic upgrade head
   ```

5. **Study Exposure**
   - New studies are hidden by default
   - Use `notebooks/06_update_qc.ipynb` to:
     - Toggle study visibility
     - Add QC comments
     - Update study metadata

### Query Interface

The database provides several categories of query functions through `queries.py`:

#### Data Access Patterns

**Metadata Queries**
- Fetch study, sample, and contrast information
- Access QC status and pipeline tracking
- Retrieve experimental design and conditions
- Query gene/transcript annotations

**Expression Data**
- Query differential expression results with flexible filtering
- Access raw and normalized expression measurements
- Support both gene and transcript-level analysis
- Calculate expression statistics across sample groups

**Analysis Support**
- Generate AnnData objects for downstream analysis
- Build summary statistics and overview tables
- Create visualization-ready data structures

See function documentation in `queries.py` for detailed parameter options and examples.


## Development Notes

### TODO
* Deprecate usage of the qc file. This is obsolete given the development of the studyqueue table. Need to automatically populate study qc from studyqueue table to study table.

