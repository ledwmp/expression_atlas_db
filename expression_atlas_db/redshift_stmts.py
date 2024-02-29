from expression_atlas_db import base

delete_stmt_all = "\n".join(
    [f"DROP TABLE if exists {t};" for t in base.Base.metadata.tables.keys()]
)
delete_stmt_sequenceregion = "\n".join(
    [
        f"DROP TABLE if exists {t};"
        for t in list(base.Base.metadata.tables.keys())[::-1]
        if t in ("sequenceregion", "gene", "transcript")
    ]
)
delete_stmt_dataset = "\n".join(
    [
        f"DROP TABLE if exists {t};"
        for t in list(base.Base.metadata.tables.keys())[::-1]
        if t not in ("sequenceregion", "gene", "transcript")
    ]
)

create_stmt_measurement = """CREATE TABLE IF NOT EXISTS differentialexpression (
	id INTEGER NOT NULL IDENTITY(1,1), 
	contrast_id INTEGER, 
	sequenceregion_id INTEGER, 
	basemean FLOAT, 
	log2foldchange FLOAT, 
	lfcse FLOAT, 
	stat FLOAT, 
	pvalue DOUBLE PRECISION, 
	log10_pvalue FLOAT,
	padj DOUBLE PRECISION, 
	log10_padj FLOAT, 
	control_mean FLOAT, 
	case_mean FLOAT, 
	PRIMARY KEY (id)
);

CREATE TABLE IF NOT EXISTS samplemeasurement (
	id INTEGER NOT NULL IDENTITY(1,1), 
	sample_id INTEGER, 
	sequenceregion_id INTEGER, 
	counts FLOAT, 
	normed_counts FLOAT, 
	tpm FLOAT, 
	normed_counts_transform FLOAT, 
	PRIMARY KEY (id)
);"""
