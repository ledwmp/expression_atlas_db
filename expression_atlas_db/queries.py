from typing import List, Union, Dict, Tuple, Callable
import pandas as pd
import numpy as np

from sqlalchemy import select, func

from expression_atlas_db import base


def unpack_fields(
    fields: Dict,
    keep: Union[List[str], Callable, None] = lambda x: x.startswith("sample_condition")
    | x.startswith("sample_type"),
) -> pd.Series:
    """ """
    if callable(keep):
        keys, values = zip(*[f for f in fields.items() if keep(f[0])])
    else:
        keys, values = zip(*[f for f in fields.items() if not keep or f[0] in keep])
    return pd.Series(values, index=keys)


def fetch_studies(session: base._Session, public: bool = True) -> pd.DataFrame:
    """ """
    query = select(base.Study)

    if public:
        query = query.filter(base.Study.public == True)

    studies = pd.read_sql(query, session.bind)

    return studies


def fetch_contrasts(
    session: base._Session, studies: Union[List[str], None] = None, public: bool = True
) -> pd.DataFrame:
    """ """
    query = select(base.Contrast, base.Study).join(
        base.Study, base.Contrast.study_id == base.Study.id
    )

    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))

    if public:
        query = query.filter(base.Study.public == True)

    contrasts_df = pd.read_sql(query, session.bind)

    return contrasts_df


def fetch_sequenceregions(
    session: base._Session,
    sequenceregions: Union[List[str], None] = None,
    sequenceregions_type: Union[str, None] = None,
) -> pd.DataFrame:

    transcript_query = select(base.Transcript)
    gene_query = select(base.Gene)

    if sequenceregions:
        transcript_query = transcript_query.filter(
            base.Transcript.transcript_id.in_(sequenceregions)
        )
        gene_query = gene_query.filter(base.Gene.gene_id.in_(sequenceregions))

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
    """ """

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
    """ """
    query = select(base.Sample, base.Study.velia_id).join(
        base.Study, base.Study.id == base.Sample.study_id
    )
    if studies:
        query = query.filter(base.Study.velia_id.in_(studies))
    if public:
        query = query.filter(base.Study.public == True)

    samples = pd.read_sql(query, session.bind)
    samples = pd.concat(
        [samples, samples["fields"].apply(lambda x: unpack_fields(x, keep_fields))],
        axis=1,
    )

    samples = samples.loc[
        :,
        ~samples.columns.str.startswith("id")
        & ~samples.columns.str.startswith("type")
        & ~samples.columns.str.startswith("study_id")
        & (samples.columns != "fields"),
    ]

    return samples


def query_differentialexpression(
    session: base._Session,
    session_redshift: base._Session,
    studies: Union[List[str], None] = None,
    contrasts: Union[List[str], None] = None,
    sequenceregions: Union[List[str], None] = None,
    sequenceregions_type: Union[str, None] = None,
    log10_padj_threshold: Union[float, None] = np.log10(0.5),
    log2_fc_threshold: Union[float, None] = np.log2(2.0),
    mean_threshold: Union[float, None] = 4.0,
    public: bool = True,
) -> pd.DataFrame:
    """ """
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
) -> pd.DataFrame:
    """ """
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
    """ """
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
