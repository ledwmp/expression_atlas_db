from typing import List, Union, Dict, Tuple, Callable, Any, Optional
import pandas as pd
import numpy as np

from sqlalchemy import select, func, or_

from expression_atlas_db import base
from expression_atlas_db.utils import MetaDataFetcher


def unpack_fields(
    fields: Dict[str, Any],
    keep: Union[List[str], Callable, None] = lambda x: x.startswith("sample_condition")
    | x.startswith("sample_type"),
) -> pd.Series:
    """Unpacks json fields column into series.

    Args:
        fields (Dict[str,Any]): dictionariy stored in the fields column as json,
            needs to be unpacked to columns.
        keep (Union[List[str], Callable, None]):
            Filter to keep specific keys in fields. Can be list of columns, or a function.
    Returns:
        (pd.Series): unpacked series with filtered index with keep.
    """
    if callable(keep):
        keys, values = zip(*[f for f in fields.items() if keep(f[0])])
    else:
        keys, values = zip(*[f for f in fields.items() if not keep or f[0] in keep])
    return pd.Series(values, index=keys)


def fetch_studies(
    session: base._Session, studies: Union[List[str], None] = None, public: bool = True
) -> pd.DataFrame:
    """Queries against the study table, returns a dataframe of the table.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        public (bool): Filter for studies without public flag set.
    Returns:
        studies_df (pd.DataFrame): studies dataframe.
    """
    query = select(base.Study)

    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))

    if public:
        query = query.filter(base.Study.public == True)

    studies_df = pd.read_sql(query, session.bind)

    return studies_df


def fetch_contrasts(
    session: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    public: bool = True,
) -> pd.DataFrame:
    """Queries against the contrast table, returns a dataframe of the table.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        contrasts (Union[List[str],None]): List of contrast_names to query against db.
        public (bool): Filter for studies without public flag set.
    Returns:
        contrasts_df (pd.DataFrame): contrasts dataframe.
    """
    query = select(base.Contrast, base.Study).join(
        base.Study, base.Contrast.study_id == base.Study.id
    )

    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))

    if contrasts:
        query = query.filter(base.Contrast.conrast_name.in_(contrasts))

    if public:
        query = query.filter(base.Study.public == True)

    contrasts_df = pd.read_sql(query, session.bind)

    return contrasts_df


def fetch_sequenceregions(
    session: base._Session,
    sequenceregions: Union[List[str], None] = None,
    sequenceregions_type: Union[str, None] = None,
    assembly_id: Union[str, None] = None,
    exact_id_match: bool = True,
) -> pd.DataFrame:
    """Queries against the sequenceregion/gene/transcript tables, returnes a dataframe of the table.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        sequenceregions (Union[List[str],None]): List of transcript_ids or gene_ids to query against database.
            These are the text ids, not the column ids from expression_atlas_db.
        sequenceregions_type (Union[str,None]): One of "transcript", "gene", or None.
            Filter query on either table or return all.
        assembly_id (Union[str, None]): VeliaDB assembly_id or public assembly id.
    Returns:
        sequenceregions_df (pd.DataFrame): sequenceregions dataframe.
    """

    transcript_query = select(base.Transcript)
    gene_query = select(base.Gene)

    if sequenceregions:
        if exact_id_match:
            transcript_query = transcript_query.filter(
                base.Transcript.transcript_id.in_(sequenceregions)
            )
            gene_query = gene_query.filter(base.Gene.gene_id.in_(sequenceregions))
        else:
            conditions_gene = [
                base.Gene.gene_id.ilike(f"{g}%") for g in sequenceregions
            ]
            conditions_transcript = [
                base.Transcript.transcript_id.ilike(f"{t}%") for t in sequenceregions
            ]
            transcript_query = transcript_query.filter(or_(*conditions_transcript))
            gene_query = gene_query.filter(or_(*conditions_gene))

    if assembly_id:
        transcript_query = transcript_query.filter(
            base.Transcript.assembly_id == assembly_id
        )
        gene_query = gene_query.filter(base.Gene.assembly_id == assembly_id)

    transcript_sequenceregions_df = pd.read_sql(transcript_query, session.bind)
    transcript_sequenceregions_df.drop("gene_id", inplace=True, axis=1)
    gene_sequenceregions_df = pd.read_sql(gene_query, session.bind)

    sequenceregions_df = pd.concat(
        [transcript_sequenceregions_df, gene_sequenceregions_df], axis=0
    )
    sequenceregions_df.set_index("id", inplace=True)

    if sequenceregions_type == "transcript":
        sequenceregions_df = sequenceregions_df[
            ~sequenceregions_df["transcript_id"].isna()
        ]
    elif sequenceregions_type == "gene":
        sequenceregions_df = sequenceregions_df[~sequenceregions_df["gene_id"].isna()]

    return sequenceregions_df


def fetch_samplecontrasts(
    session: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    keep_fields: Union[List[str], Callable, None] = lambda x: x.startswith(
        "sample_condition"
    )
    | x.startswith("sample_type"),
    public: bool = True,
) -> pd.DataFrame:
    """Queries against the samplecontrast table, returns a dataframe of contrasts merged
    with the samples used to create the contrast.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        contrasts (Union[List[str],None]): List of contrast_names to query against db.
        keep_fields (Union[List[str], Callable, None]): Passed directly to unpack_fields.
            Columns to keep or callable to filter columns.
    Returns:
        samplecontrast_df (pd.DataFrame): samplecontrasts dataframe.
    """

    query = (
        select(
            base.Contrast,
            base.Study.velia_id,
            base.SampleContrast.contrast_side,
            base.Sample,
        )
        .join(base.Study, base.Study.id == base.Contrast.study_id)
        .join(base.SampleContrast, base.Contrast.id == base.SampleContrast.contrast_id)
        .join(base.Sample, base.SampleContrast.sample_id == base.Sample.id)
    )
    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))
    if public:
        query = query.filter(base.Study.public == True)
    if contrasts:
        query = query.filter(base.Contrast.contrast_name.in_(contrasts))

    samplecontrasts_df = pd.read_sql(query, session.bind)

    samplecontrasts_df = pd.concat(
        [
            samplecontrasts_df,
            samplecontrasts_df["fields"].apply(lambda x: unpack_fields(x, keep_fields)),
        ],
        axis=1,
    )

    samplecontrasts_df = samplecontrasts_df.loc[
        :,
        ~samplecontrasts_df.columns.str.startswith("id")
        & ~samplecontrasts_df.columns.str.startswith("type")
        & ~samplecontrasts_df.columns.str.startswith("study_id")
        & (samplecontrasts_df.columns != "fields"),
    ]

    return samplecontrasts_df


def fetch_samples(
    session: base._Session,
    studies: Union[List[str], None] = None,
    keep_fields: Union[List[str], Callable, None] = lambda x: x.startswith(
        "sample_condition"
    ),
    public: bool = True,
) -> pd.DataFrame:
    """Queries against the samples table, returns a dataframe of the samples.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        keep_fields (Union[List[str], Callable, None]): Passed directly to unpack_fields.
            Columns to keep or callable to filter columns.
    Returns:
        samples_df (pd.DataFrame): samples dataframe.

    """
    query = select(base.Sample, base.Study.velia_id).join(
        base.Study, base.Study.id == base.Sample.study_id
    )
    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))

    if public:
        query = query.filter(base.Study.public == True)

    samples_df = pd.read_sql(query, session.bind)
    samples_df = pd.concat(
        [
            samples_df,
            samples_df["fields"].apply(lambda x: unpack_fields(x, keep_fields)),
        ],
        axis=1,
    )

    samples_df = samples_df.loc[
        :,
        ~samples_df.columns.str.startswith("id")
        & ~samples_df.columns.str.startswith("type")
        & ~samples_df.columns.str.startswith("study_id")
        & (samples_df.columns != "fields"),
    ]

    return samples_df


def fetch_studyqueue(
    session: base._Session,
    public: bool = True,
) -> pd.DataFrame:
    """Queries against the studyqueue table, returns a dataframe of the studies in queue.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        public (bool): Filter for studies without public flag set.
    Returns:
        studyqueues_df (pd.DataFrame): samples dataframe.
    """
    query = select(
        base.StudyQueue,
        base.Study.id.label("study_id"),
        base.Study.public.label("data_released"),
    ).join(
        base.Study,
        base.StudyQueue.study_id == base.Study.id,
        isouter=True,
    )
    if public:
        query = query.filter(base.StudyQueue.public == True)

    return pd.read_sql(
        query,
        con=session.bind,
    )


def update_studyqueue(
    session: base._Session,
    update_rows: pd.DataFrame,
    update_columns: List[str] = [
        "request",
        "status",
        "technology",
        "quality",
    ],
) -> pd.DataFrame:
    """Update studyqueue table with metadata from editable dataframe.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        update_rows (pd.DataFrame): Dataframe of rows to be updated in studyqueue.
        update_columns (List[str]): List of column names to restrict update to.
    Returns:
        (pd.DataFrame): Updated rows in dataframe.
    """
    for i, r in update_rows.iterrows():
        sq = (
            session.query(base.StudyQueue).filter(base.StudyQueue.id == r["id"]).first()
        )
        for c in update_columns:
            setattr(sq, c, r.get(c))
        if r.get("delete", False) == True:
            setattr(sq, "public", False)
    session.commit()
    return pd.read_sql(
        select(base.StudyQueue).filter(
            base.StudyQueue.id.in_(update_rows["id"].tolist())
        ),
        con=session.bind,
    )


def submit_studyqueue(
    session: base._Session,
    srp_id: str,
    request: str,
    technology: str,
    geo_id: Union[str, None] = None,
) -> Tuple[bool, Optional[pd.DataFrame], Optional[str]]:
    """Submit a study to the queue table.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        srp_id (str): SRP project id.
        request (str): Reason for requesting dataset.
        technology (str): Technology type for queue. Must be one of {'BULK', '10X-3', '10X-5', 'SMART-SEQ'}.
        geo_id (str): GEO ID for project.

    Returns:
        (Tuple[bool, Optional[pd.DataFrame], Optional[str]]): Study exists, dataframe with results, message string.

    """
    if technology not in (
        "BULK",
        "10X-3",
        "10X-5",
        "SMART-SEQ",
    ):
        return (
            False,
            f"Technology must be one of {'BULK', '10X-3', '10X-5', 'SMART-SEQ'}, is: {srp_id}.",
        )
    studyqueue = (
        session.query(base.StudyQueue)
        .filter((base.StudyQueue.srp_id == srp_id) | (base.StudyQueue.geo_id == geo_id))
        .first()
    )

    if not studyqueue:
        meta = MetaDataFetcher(srp_id, [])
        if (meta.srp_id and geo_id and meta.geo_id == geo_id) or (
            not geo_id and meta.srp_id
        ):
            studyqueue = base.StudyQueue(
                velia_id=srp_id,
                geo_id=geo_id,
                srp_id=srp_id,
                pmid=meta.pmids,
                title=meta.project_title,
                description=meta.project_summary,
                public=True,
                processed=False,
                status="QUEUED",
                technology=technology,
                request=request,
            )
            session.add(studyqueue)
            session.commit()
            return (
                False,
                pd.read_sql(
                    select(base.StudyQueue).filter(base.StudyQueue.id == studyqueue.id),
                    con=session.bind,
                ),
            )
        elif geo_id and meta.geo_id != geo_id:
            return (
                False,
                f"Submitted geo id and geo id associated with {srp_id} don't match.",
            )
        else:
            return (False, f"Cannot request metadata for {srp_id}.")

    else:
        return (
            True,
            pd.read_sql(
                select(base.StudyQueue).filter(base.StudyQueue.id == studyqueue.id),
                con=session.bind,
            ),
        )


def query_differentialexpression(
    session: base._Session,
    session_redshift: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    sequenceregions: Union[List[str], None] = None,
    sequenceregions_type: Union[str, None] = None,
    log10_padj_threshold: Union[float, None] = np.log10(0.05),
    log2_fc_threshold: Union[float, None] = np.log2(2.0),
    mean_threshold: Union[float, None] = 4.0,
    public: bool = True,
    exact_id_match: bool = True,
) -> pd.DataFrame:
    """Queries against both databases to fetch entries out of the differentialexpression
    table. First, queries the contrast and sequenceregion tables in the postgres/sqlite db,
    then uses the contrast_ids and sample_ids to query the differentialexpression table in
    redshift.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        contrasts (Union[List[str],None]): List of contrast_names to query against db.
        sequenceregions (Union[List[str],None]): List of transcript_ids or gene_ids to query against database.
            These are the text ids, not the column ids from expression_atlas_db.
        sequenceregions_type (Union[str,None]): One of "transcript", "gene", or None.
            Filter query on either table or return all.
        log10_padj_threshold (Union[float,None]): Keep rows log10_padj < threshold.
        log2_fc_threshold (Union[float,None]): Keep rows abs(log2(fc)) > threshold.
        mean_threshold (Union[float,None]): Keep rows mean(normed_transformed_count) in case or control.
        public (bool): Filter for studies without public flag set.
        exact_id_match (bool): Whether to allow matches to sequenceregions with prefixes.
    Returns:
        differentialexpression_df (pd.DataFrame): differentialexpression dataframe.
    """
    studies_query = select(base.Contrast, base.Study.velia_id).join(
        base.Study, base.Contrast.study_id == base.Study.id
    )

    if studies:
        studies_query = studies_query.filter(base.Study.velia_id.in_(studies))

    if public:
        studies_query = studies_query.filter(base.Study.public == True)

    if contrasts:
        studies_query = studies_query.filter(base.Contrast.contrast_name.in_(contrasts))

    studies_df = pd.read_sql(studies_query, session.bind)
    studies_df.set_index("id", inplace=True)

    sequenceregions_df = fetch_sequenceregions(
        session,
        sequenceregions=sequenceregions,
        sequenceregions_type=sequenceregions_type,
        exact_id_match=exact_id_match,
    )

    differentialexpression_query = select(base.DifferentialExpression).filter(
        base.DifferentialExpression.contrast_id.in_(studies_df.index)
    )

    if sequenceregions:
        differentialexpression_query = differentialexpression_query.filter(
            base.DifferentialExpression.sequenceregion_id.in_(sequenceregions_df.index)
        )

    if log10_padj_threshold:
        differentialexpression_query = differentialexpression_query.filter(
            base.DifferentialExpression.log10_padj <= log10_padj_threshold
        )

    if mean_threshold:
        differentialexpression_query = differentialexpression_query.filter(
            (base.DifferentialExpression.case_mean >= mean_threshold)
            | (base.DifferentialExpression.control_mean >= mean_threshold)
        )

    if log2_fc_threshold:
        differentialexpression_query = differentialexpression_query.filter(
            func.abs(base.DifferentialExpression.log2foldchange) >= log2_fc_threshold
        )

    differentialexpression_df = pd.read_sql(
        differentialexpression_query, session_redshift.bind
    )
    differentialexpression_df = differentialexpression_df.merge(
        studies_df[["velia_id", "contrast_name"]],
        left_on="contrast_id",
        right_index=True,
        how="left",
    ).merge(
        sequenceregions_df[["transcript_id", "gene_id", "type"]],
        left_on="sequenceregion_id",
        right_index=True,
    )

    differentialexpression_df = differentialexpression_df.loc[
        :,
        ~differentialexpression_df.columns.str.startswith("id")
        & ~differentialexpression_df.columns.str.startswith("contrast_id")
        & ~differentialexpression_df.columns.str.startswith("sequenceregion_id")
        & ~differentialexpression_df.columns.str.startswith("type"),
    ]

    return differentialexpression_df


def query_samplemeasurement(
    session: base._Session,
    session_redshift: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    samples: Union[List[str], None] = None,
    sequenceregions: Union[List[str], None] = None,
    sequenceregions_type: Union[str, None] = None,
    public: bool = True,
    exact_id_match: bool = True,
) -> pd.DataFrame:
    """Queries against both databases to fetch entries out of the differentialexpression
    table. First, queries the contrast and sequenceregion tables in the postgres/sqlite db,
    then uses the contrast_ids and sample_ids to query the differentialexpression table in
    redshift.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        contrasts (Union[List[str],None]): List of contrast_names to query against db.
        sequenceregions (Union[List[str],None]): List of transcript_ids or gene_ids to query against database.
            These are the text ids, not the column ids from expression_atlas_db.
        sequenceregions_type (Union[str,None]): One of "transcript", "gene", or None.
            Filter query on either table or return all.
        exact_id_match (bool): Whether to allow matches to sequenceregions with prefixes.
    Returns:
        samplemeasurement_df (pd.DataFrame): samplemeasurements dataframe.
    """
    samples_query = select(base.Sample, base.Study.velia_id).join(
        base.Sample, base.Sample.study_id == base.Study.id
    )
    if studies:
        samples_query = samples_query.filter(base.Study.velia_id.in_(studies))

    if public:
        samples_query = samples_query.filter(base.Study.public == True)

    if samples:
        samples_query = samples_query.filter(base.Sample.srx_id.in_(samples))

    if contrasts:
        contrast_subquery = (
            select(
                base.Study.velia_id,
                base.Contrast.contrast_name,
                base.SampleContrast.contrast_side,
                base.SampleContrast.sample_id,
            )
            .filter(base.Study.velia_id.in_(studies))
            .filter(base.Contrast.contrast_name.in_(contrasts))
            .join(base.Contrast, base.Contrast.study_id == base.Study.id)
            .join(
                base.SampleContrast, base.Contrast.id == base.SampleContrast.contrast_id
            )
            .subquery()
        )
        samples_subquery = samples_query.subquery()
        samples_query = select(contrast_subquery, samples_subquery).join(
            contrast_subquery, samples_subquery.c.id == contrast_subquery.c.sample_id
        )

    samples_df = pd.read_sql(samples_query, session.bind)
    samples_df.set_index("id", inplace=True)

    sequenceregions_df = fetch_sequenceregions(
        session,
        sequenceregions=sequenceregions,
        sequenceregions_type=sequenceregions_type,
        exact_id_match=exact_id_match,
    )

    samplemeasurement_query = select(base.SampleMeasurement).filter(
        base.SampleMeasurement.sample_id.in_(samples_df.index)
    )

    if sequenceregions:
        samplemeasurement_query = samplemeasurement_query.filter(
            base.SampleMeasurement.sequenceregion_id.in_(sequenceregions_df.index)
        )

    samplemeasurement_df = pd.read_sql(samplemeasurement_query, session_redshift.bind)
    samplemeasurement_df = samplemeasurement_df.merge(
        samples_df[
            ["velia_id", "srx_id", "atlas_group", "fields"]
            + ([] if not contrasts else ["contrast_name", "contrast_side"])
        ],
        left_on="sample_id",
        right_index=True,
        how="left",
    ).merge(
        sequenceregions_df[["transcript_id", "gene_id", "type"]],
        left_on="sequenceregion_id",
        right_index=True,
    )

    samplemeasurement_df = samplemeasurement_df.loc[
        :,
        ~samplemeasurement_df.columns.str.startswith("id")
        & ~samplemeasurement_df.columns.str.startswith("type")
        & ~samplemeasurement_df.columns.str.startswith("study_id")
        & ~samplemeasurement_df.columns.str.startswith("sample_id"),
    ]

    samplemeasurement_df = pd.concat(
        [
            samplemeasurement_df,
            samplemeasurement_df["fields"].apply(
                lambda x: unpack_fields(
                    x,
                    lambda x: x
                    in (
                        "sample_condition_1",
                        "sample_condition_2",
                        "sample_type_1",
                        "sample_type_2",
                    ),
                )
            ),
        ],
        axis=1,
    )

    return samplemeasurement_df


def build_contrast_metatable(
    session: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    public: bool = True,
) -> pd.DataFrame:
    """Helper function to assemble the dashboard contrast metatable. Calls fetch_samplecontrasts
    and coerces table to format required by dashboard.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        studies (Union[List[str],None]): List of velia_id studies to query against db.
        contrasts (Union[List[str],None]): List of contrast_names to query against db.
        public (bool): Filter for studies without public flag set.
    Returns:
        samplecontrast_df (pd.DataFrame): dashboard contrast metadata dataframe.
    """
    samplecontrast_df = fetch_samplecontrasts(
        session,
        studies=studies,
        contrasts=contrasts,
        public=public,
    ).fillna("")

    samplecontrast_df.drop(
        [
            "sample_condition_key",
            "left_condition_level",
            "right_condition_level",
            "left_condition_display",
            "right_condition_display",
            "atlas_group",
        ],
        axis=1,
        inplace=True,
    )

    samplecontrast_df = samplecontrast_df.astype(str)

    samplecontrast_df.insert(
        1,
        "sample_n",
        samplecontrast_df.groupby(["velia_id", "contrast_name", "contrast_side"])[
            "srx_id"
        ]
        .transform(len)
        .astype(str),
    )
    samplecontrast_df.drop("srx_id", inplace=True, axis=1)
    samplecontrast_df = samplecontrast_df.groupby(
        ["velia_id", "contrast_name", "contrast_side"],
    ).agg(lambda x: ",".join(x.unique()))
    for c in samplecontrast_df.columns[
        samplecontrast_df.columns.str.startswith("sample")
    ]:
        samplecontrast_df[f"left_{c}"] = samplecontrast_df.loc[
            samplecontrast_df.index.get_level_values("contrast_side") == "left"
        ][c]
        samplecontrast_df[f"right_{c}"] = samplecontrast_df.loc[
            samplecontrast_df.index.get_level_values("contrast_side") == "right"
        ][c]

    samplecontrast_df.reset_index(inplace=True, drop=False)
    samplecontrast_df.drop(
        samplecontrast_df.columns[
            samplecontrast_df.columns.str.startswith("sample")
            | samplecontrast_df.columns.str.startswith("created_at")
            | samplecontrast_df.columns.str.startswith("alembic_id")
        ].tolist()
        + ["contrast_side"],
        inplace=True,
        axis=1,
    )
    samplecontrast_df.index = samplecontrast_df.apply(
        lambda x: f"{x['velia_id']} -- {x['contrast_name']}", axis=1
    )
    samplecontrast_df = (
        samplecontrast_df.groupby(
            [samplecontrast_df.index, "velia_id", "contrast_name"]
        )
        .agg(sum)
        .reset_index(["velia_id", "contrast_name"], drop=False)
    )
    return samplecontrast_df


def query_percentile_group(
    session: base._Session,
    session_redshift: base._Session,
    sample_ids: List[int],
    percentile_levels: List[float] = [0.5, 0.75, 0.90],
    sequenceregion_ids: Union[List[int], None] = None,
    aggregate_column: str = "tpm",
) -> pd.DataFrame:
    """
    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        session_redshift (base._Session): SQLAlchemy session object to redshift db.
        sample_ids (Union[List[str],None]): List of sample_ids to query against db.
        percentile_levels (List[float]): List of percentile levels to query.
        sequenceregion_ids (Union[List[int],None]): List of optional sequenceregion_ids to query.
        aggregate_column (str): Column to aggregate percentiles over.
    Returns:
        percentile_df (pd.DataFrame): Dataframe with sequenceregion_ids and aggregated expression percentiles.
    """

    if aggregate_column not in vars(base.SampleMeasurement).keys():
        raise KeyError("Aggregate_column does not exist in samplemeasurement table.")

    query = select(
        base.SampleMeasurement.sequenceregion_id,
        *[
            func.percentile_disc(i)
            .within_group(vars(base.SampleMeasurement)[aggregate_column])
            .over(partition_by=base.SampleMeasurement.sequenceregion_id)
            .label(f"perc_{i}")
            for i in percentile_levels
        ],
    ).filter(base.SampleMeasurement.sample_id.in_(sample_ids))

    if sequenceregion_ids:
        query = query.filter(
            base.SampleMeasurement.sequenceregion_id.in_(sequenceregion_ids)
        )

    percentile_df = pd.read_sql(query.distinct(), session_redshift.bind)

    return percentile_df


def build_expression_atlas_summary_dfs(
    session: base._Session,
    public: bool = True,
    ambiguous_tisues: List[str] = ["COLON", "BLOOD", "BONE_MARROW", "BRAIN"],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Builds summary study, contrast, and sample tables from expression atlas.

    Args:
        session (base._Session): SQLAlchemy session object to the main postgres/sqlite db.
        public (bool): Filter for studies without public flag set.
        ambiguous_tissues (List[str]): List of tissues (sample_type_1) where we want to transfer higher-granularity labels.
    Returns:
        studies_df, contrasts_df, samples_df (Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]): dfs.
    """
    # Pull studies df.
    studies_df = fetch_studies(session, public=public)

    # Pull and munge metadata samples df.
    samples_df = fetch_samples(
        session,
        keep_fields=lambda x: x.startswith(
            (
                "sample_type_1",
                "sample_type_2",
                "sample_condition_1",
                "sample_condition_2",
                "disease_",
            )
        ),
        public=public,
    )
    samples_df["primary_disease"] = samples_df.apply(
        lambda x: [
            c
            for c in x[x.index.str.startswith("disease_")]
            if not pd.isna(c) and c != "NA"
        ],
        axis=1,
    )
    samples_df["primary_tissue"] = samples_df.apply(
        lambda x: (
            x["sample_type_2"].upper()
            if not pd.isna(x["sample_type_2"])
            and x["sample_type_1"] in ambiguous_tisues
            else x["sample_type_1"].upper()
        ),
        axis=1,
    )
    samples_df["primary_condition"] = samples_df["sample_condition_1"].copy()

    # Pull samplecontrasts link table.
    contrasts_df = fetch_samplecontrasts(session, public=public)

    # Merge metadata into studies, contrasts.
    studies_df = (
        studies_df.loc[:, ["velia_id", "pmid", "title", "description"]]
        .merge(
            samples_df.loc[:, ["velia_id", "primary_tissue", "primary_disease"]],
            on="velia_id",
            how="left",
        )
        .explode("primary_disease")
        .drop_duplicates()
        .reset_index(drop=True)
    )

    contrasts_df = (
        contrasts_df.loc[
            contrasts_df["contrast_side"] == "left",
            ["velia_id", "contrast_name", "srx_id"],
        ]
        .merge(
            samples_df.loc[
                :, ["velia_id", "primary_tissue", "primary_disease", "srx_id"]
            ],
            on=["velia_id", "srx_id"],
            how="left",
        )
        .loc[:, ["velia_id", "contrast_name", "primary_tissue", "primary_disease"]]
        .explode("primary_disease")
        .drop_duplicates()
        .reset_index(drop=True)
    )

    samples_df = samples_df.loc[
        :,
        [
            "srx_id",
            "velia_id",
            "primary_tissue",
            "primary_disease",
            "primary_condition",
        ],
    ]

    return studies_df, contrasts_df, samples_df
