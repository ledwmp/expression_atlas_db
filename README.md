# Expression Atlas Database

### TODO:
* Swap GTF feature load to veliadb. v0c/v1 veliadb ids broken. Will be better to migrate to load from veliadb itself.
* Deprecate usage of of the qc file. This is obsolete given the development of the studyqueue table. Need to automatically populate study qc from studyqueue table to study table.

### Steps to trigger study ingestion:
* Not best use case for alembic, but all database loads are handled via an alembic revision. Overkill, but makes managing the loads easy and can wipe if something goes wrong.  
1. Copy study directory into s3 at: `s3://velia-piperuns-dev/expression_atlas/v1/`. The study (ex. SRP000001) directory should have structure like:
    ```
    SRP000001  
    |   ...
    └───fetchngs_output
    │   │   ...
    └───rnaseq_output
        │   ...
    └───de_results
        │   SRP000001_dds_transcript.h5_ad
        │   SRP000001_dds_gene.h5_ad
        |   ...
    ``` 
    The de_results folder and anndatas are the only required files for ingestion. There are a few strict requirements for the folder/anndatas:
    1. obs must be indexed on the srx accession id.
    2. var must be indexed on the gene_id or transcript_id attributes from the veliadb gtf.
    3. Layers exist: `["counts", "normed_counts", "raw_tpm", "normed_counts_transform"]`
    4. uns must have dictionary with name "contrasts" that is keyed on a contrast name with values an iterable containing: `[<factor name>, <test factor>, <reference factor>]`, where factor name is a column name in obs, and the factors are categorical levels of that column.
    5. uns must have a dictionary with name "stat_results" that is keyed on contrast name from 4. with values as a dataframe with structure:
        1. indexed with similar index as var from 2.
        2. columns: `["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "-log10_padj"]` and an additional two columns with case and control means of normed_counts_transform named like: `<factor name>_<factor level>_meannormedcounts`.
    6. Folder/project should be named (ex. SRP000001) as an SRA/ENA project id. Metadata are automatically fetched from ENA, so this is needed for that step.
    7. The anndatas must have the exact suffix described above.

    All of the above should be automatically formatted in the `de_processing` notebooks 01-04 of the expression_atlas repo.  
2. At base level of project directory, begin database load with below code block.
    ```
    alembic revision -m 'add'
    alembic upgrade head
    ```
    This walks over the expression atlas piperuns bucket and loads studies. This is simply calling load_db.add_studies followed by load_db.write_studies_qc while tagging the added rows in the dataset table with an update id. One could easily call these functions out of a notebook. 
3. By default, studies are hidden until exposed manually. Run notebook `notebooks/06_update_qc.ipynb` and toggle public to true for datsets that need to be exposed to dashboard. Additionally, qc comments can be manually added here or merged in from queue table. Finally, this notebook calls load_db.update_studies_qc to transfer this updated qc table to studies table.

