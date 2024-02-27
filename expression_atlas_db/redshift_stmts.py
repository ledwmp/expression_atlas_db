from expression_atlas_db import base

delete_stmt_all = '\n'.join([f'DROP TABLE if exists {t};' for t in base.Base.metadata.tables.keys()])
delete_stmt_sequenceregion = '\n'.join([f'DROP TABLE if exists {t};' for t in list(base.Base.metadata.tables.keys())[::-1] \
                                                    if t in ('sequenceregion','gene','transcript')])
delete_stmt_dataset = '\n'.join([f'DROP TABLE if exists {t};' for t in list(base.Base.metadata.tables.keys())[::-1] \
                                                    if t not in ('sequenceregion','gene','transcript')])

# create_stmt_dataset = """CREATE TABLE IF NOT EXISTS dataset (
# 	id INTEGER NOT NULL IDENTITY(1,1), 
# 	type VARCHAR(30), 
# 	PRIMARY KEY (id)
# );

# CREATE TABLE IF NOT EXISTS study (
# 	id INTEGER NOT NULL, 
# 	velia_id VARCHAR(20) NOT NULL, 
# 	geo_id VARCHAR(20), 
# 	srp_id VARCHAR(20), 
# 	bio_id VARCHAR(20), 
# 	pmid INTEGER, 
# 	is_public BOOLEAN, 
# 	quality VARCHAR(20), 
# 	timestamps VARCHAR(100), 
# 	sizes VARCHAR(100), 
# 	title VARCHAR(200), 
# 	description TEXT, 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES dataset (id)
# );

# CREATE TABLE IF NOT EXISTS sample (
# 	id INTEGER NOT NULL, 
# 	study_id INTEGER, 
# 	srx_id VARCHAR(20) NOT NULL, 
# 	atlas_group VARCHAR(200) NOT NULL, 
# 	fields SUPER, 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES dataset (id), 
# 	FOREIGN KEY(study_id) REFERENCES study (id)
# );

# CREATE TABLE IF NOT EXISTS contrast (
# 	id INTEGER NOT NULL, 
# 	study_id INTEGER, 
# 	contrast_name VARCHAR(200) NOT NULL, 
# 	sample_condition_key VARCHAR(200), 
# 	left_condition_level VARCHAR(200), 
# 	right_condition_level VARCHAR(200), 
# 	left_condition_display VARCHAR(200), 
# 	right_condition_display VARCHAR(200), 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES dataset (id), 
# 	FOREIGN KEY(study_id) REFERENCES study (id)
# );

# CREATE TABLE IF NOT EXISTS samplecontrast (
# 	id INTEGER NOT NULL, 
# 	sample_id INTEGER, 
# 	contrast_id INTEGER, 
# 	study_id INTEGER, 
# 	sample_condition_key VARCHAR(200), 
# 	condition_level VARCHAR(200), 
# 	condition_display VARCHAR(200), 
# 	contrast_side VARCHAR(10), 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES dataset (id), 
# 	FOREIGN KEY(sample_id) REFERENCES sample (id), 
# 	FOREIGN KEY(contrast_id) REFERENCES contrast (id), 
# 	FOREIGN KEY(study_id) REFERENCES study (id)
# );"""

# create_stmt_measurement = """CREATE TABLE IF NOT EXISTS differentialexpression (
# 	id INTEGER NOT NULL IDENTITY(1,1), 
# 	contrast_id INTEGER, 
# 	sequenceregion_id INTEGER, 
# 	basemean FLOAT, 
# 	log2foldchange FLOAT, 
# 	lfcse FLOAT, 
# 	stat FLOAT, 
# 	pvalue FLOAT, 
# 	padj FLOAT, 
# 	log10_padj FLOAT, 
# 	control_mean FLOAT, 
# 	case_mean FLOAT, 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(contrast_id) REFERENCES contrast (id), 
# 	FOREIGN KEY(sequenceregion_id) REFERENCES sequenceregion (id)
# );

# CREATE TABLE IF NOT EXISTS samplemeasurement (
# 	id INTEGER NOT NULL IDENTITY(1,1), 
# 	sample_id INTEGER, 
# 	sequenceregion_id INTEGER, 
# 	counts FLOAT, 
# 	normed_counts FLOAT, 
# 	tpm FLOAT, 
# 	normed_counts_transform FLOAT, 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(sample_id) REFERENCES sample (id), 
# 	FOREIGN KEY(sequenceregion_id) REFERENCES sequenceregion (id)
# );"""


create_stmt_measurement = """CREATE TABLE IF NOT EXISTS differentialexpression (
	id INTEGER NOT NULL IDENTITY(1,1), 
	contrast_id INTEGER, 
	sequenceregion_id INTEGER, 
	basemean FLOAT, 
	log2foldchange FLOAT, 
	lfcse FLOAT, 
	stat FLOAT, 
	pvalue FLOAT, 
	padj FLOAT, 
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

# create_stmt_sequenceregion = """CREATE TABLE IF NOT EXISTS sequenceregion (
# 	id INTEGER NOT NULL IDENTITY(1,1), 
# 	veliadb_id INTEGER NOT NULL, 
# 	assembly_id VARCHAR(200), 
# 	type VARCHAR(30), 
# 	PRIMARY KEY (id)
# );

# CREATE TABLE IF NOT EXISTS gene (
# 	id INTEGER NOT NULL, 
# 	gene_id VARCHAR(50), 
# 	gene_biotype VARCHAR(200), 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES sequenceregion (id)
# );

# CREATE TABLE IF NOT EXISTS transcript (
# 	id INTEGER NOT NULL, 
# 	transcript_id VARCHAR(50), 
# 	gene_id INTEGER, 
# 	gene_biotype VARCHAR(200), 
# 	PRIMARY KEY (id), 
# 	FOREIGN KEY(id) REFERENCES sequenceregion (id), 
# 	FOREIGN KEY(gene_id) REFERENCES gene (id)
# );"""