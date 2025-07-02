#!/usr/bin/env python

__author__ = "Cote Falaguera (mjfalagueramata@gmail.com)"
__date__ = "02 Jul 2025"

"""
timeseries.py: Assess the evolution over time of evidence supporting target-disease associations in the Open Targets Platform.

Useful GitHub links:
- https://github.com/opentargets/timeseries
- https://github.com/opentargets/issues/issues/2739
"""

# Setup Google Cloud machine and sumbit job:
# gcloud dataproc clusters create cf-timeseries --image-version 2.2 --region europe-west1 --master-machine-type n1-standard-2 --secondary-worker-type spot --worker-machine-type n1-standard-4 --worker-boot-disk-size 500 --autoscaling-policy=otg-etl --optional-components=JUPYTER --enable-component-gateway --project open-targets-eu-dev
# gcloud dataproc jobs submit pyspark timeseries.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"


import datetime
import os
import time
from datetime import timedelta

from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window
from pyspark.sql import types as T


# Required data
first_year = 2000
last_year = datetime.date.today().year
novelty_scale = 2  # 2 for long tailed slow decays
novelty_shift = 2  # 3 for long tailed slow decays
novelty_window = 10
max_score = 1.64  # np.sum(1 / np.arange(1, 10000 + 1) ** 2)

data_source = [
    {
        "id": "gwas_credible_sets",
        "sectionId": "gwasCredibleSets",
        "label": "GWAS associations",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,  # needs to be a float
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gwas-associations",
    },
    {
        "id": "eva",
        "sectionId": "eva",
        "label": "ClinVar",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#clinvar",
    },
    {
        "id": "gene_burden",
        "sectionId": "geneBurden",
        "label": "Gene Burden",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gene-burden",
    },
    {
        "id": "genomics_england",
        "sectionId": "genomicsEngland",
        "label": "GEL PanelApp",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#genomics-england-panelapp",
    },
    {
        "id": "gene2phenotype",
        "sectionId": "gene2Phenotype",
        "label": "Gene2phenotype",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gene2phenotype",
    },
    {
        "id": "uniprot_literature",
        "sectionId": "uniprotLiterature",
        "label": "UniProt literature",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#uniprot-literature",
    },
    {
        "id": "uniprot_variants",
        "sectionId": "uniprotVariants",
        "label": "UniProt curated variants",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#uniprot-variants",
    },
    {
        "id": "orphanet",
        "sectionId": "orphanet",
        "label": "Orphanet",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#orphanet",
    },
    {
        "id": "clingen",
        "sectionId": "clinGen",
        "label": "Clingen",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#clingen",
    },
    {
        "id": "cancer_gene_census",
        "sectionId": "cancerGeneCensus",
        "label": "Cancer Gene Census",
        "aggregation": "Somatic mutations",
        "aggregationId": "somatic_mutation",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#cancer-gene-census",
    },
    {
        "id": "intogen",
        "sectionId": "intOgen",
        "label": "IntOGen",
        "aggregation": "Somatic mutations",
        "aggregationId": "somatic_mutation",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#intogen",
    },
    {
        "id": "eva_somatic",
        "sectionId": "evaSomatic",
        "label": "ClinVar (somatic)",
        "aggregation": "Somatic mutations",
        "aggregationId": "somatic_mutation",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#clinvar-somatic",
    },
    {
        "id": "cancer_biomarkers",
        "sectionId": "cancerBiomarkers",
        "label": "Cancer Biomarkers",
        "aggregation": "Somatic mutations",
        "aggregationId": "somatic_mutation",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#cancer-biomarkers",
    },
    {
        "id": "chembl",
        "sectionId": "chembl",
        "label": "ChEMBL",
        "aggregation": "Known drug",
        "aggregationId": "known_drug",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#chembl",
    },
    {
        "id": "crispr_screen",
        "sectionId": "crispr_screen",
        "label": "CRISPR Screens",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#project-score",
    },
    {
        "id": "crispr",
        "sectionId": "crispr",
        "label": "Project Score",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#project-score",
    },
    {
        "id": "slapenrich",
        "sectionId": "slapEnrich",
        "label": "SLAPenrich",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 0.5,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#slapenrich",
    },
    {
        "id": "progeny",
        "sectionId": "progeny",
        "label": "PROGENy",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 0.5,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#slapenrich",
    },
    {
        "id": "reactome",
        "sectionId": "reactome",
        "label": "Reactome",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#reactome",
    },
    {
        "id": "sysbio",
        "sectionId": "sysBio",
        "label": "Gene signatures",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 0.5,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gene-signatures",
    },
    {
        "id": "europepmc",
        "sectionId": "europePmc",
        "label": "Europe PMC",
        "aggregation": "Literature",
        "aggregationId": "literature",
        "weight": 0.2,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#europe-pmc",
    },
    {
        "id": "expression_atlas",
        "sectionId": "expression",
        "label": "Expression Atlas",
        "aggregation": "RNA expression",
        "aggregationId": "rna_expression",
        "weight": 0.2,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#expression-atlas",
    },
    {
        "id": "impc",
        "sectionId": "impc",
        "label": "IMPC",
        "aggregation": "Animal model",
        "aggregationId": "animal_model",
        "weight": 0.2,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#impc",
    },
    # {
    #     "id": "ot_crispr",
    #     "sectionId": "otCrispr",
    #     "label": "OT CRISPR",
    #     "aggregation": "Partner-only",
    #     "aggregationId": "partner_only",
    #     "weight": 0.5,
    #     "isPrivate": True,
    #     "docsLink": "https://partner-platform.opentargets.org/projects",
    # },
    # {
    #     "id": "encore",
    #     "sectionId": "encore",
    #     "label": "ENCORE",
    #     "aggregation": "Partner-only",
    #     "aggregationId": "partner_only",
    #     "weight": 0.5,
    #     "isPrivate": True,
    #     "docsLink": "https://partner-platform.opentargets.org/projects",
    # },
    # {
    #     "id": "ot_crispr_validation",
    #     "sectionId": "validationlab",
    #     "label": "OT Validation",
    #     "aggregation": "Partner-only",
    #     "aggregationId": "partner_only",
    #     "weight": 0.5,
    #     "isPrivate": True,
    #     "docsLink": "https://partner-platform.opentargets.org/projects",
    # },
]


# Paths
ot_platform_version = "25.03"

data_path = "gs://open-targets-data-releases/{}/output/".format(ot_platform_version)
disease_file = data_path + "disease"

results_path = "gs://ot-team/cfalaguera/{}/".format(ot_platform_version)
evidence_dated_file = results_path + "/evidence_dated/"
evidence_dated_indirect_file = results_path + "/evidence_dated_indirect/"
association_by_datasource_dated_indirect_file = (
    results_path + "association_by_datasource_dated_indirect"
)
association_by_overall_dated_indirect_file = (
    results_path + "association_by_overall_dated_indirect"
)
association_by_datasource_dated_file = results_path + "association_by_datasource_dated"
association_by_overall_dated_file = results_path + "association_by_overall_dated"

# novelty_by_datasource_dated_file = results_path + "novelty_by_datasource_dated"
# novelty_by_datasource_dated_indirect_file = (
#     results_path + "novelty_by_datasource_dated_indirect"
# )

# Establish spark connection
spark = SparkSession.builder.getOrCreate()


def get_weight_for_datasource():
    """
    Returns list of data sources' weights for overall score.
    """

    weights = [[datasource["id"], datasource["weight"]] for datasource in data_source]

    return weights


# Novelty assessment functions
def get_indirect_evidence():
    """
    Propagate dated evidence accross disease ontology.
    """

    if os.path.exists(evidence_dated_indirect_file):
        pass

    else:

        (
            spark.read.parquet(evidence_dated_file)
            .join(
                spark.read.parquet(disease_file).select(
                    F.col("id").alias("diseaseId"),
                    F.explode(
                        F.array_union(F.array(F.col("id")), F.col("ancestors"))
                    ).alias(  # add descendant to ancestors list to keep it for later
                        "specificDiseaseId"
                    ),
                ),
                "diseaseId",
                "inner",
            )  # add ancestors
            .drop("diseaseId")  # drop descendants
            .withColumnRenamed("specificDiseaseId", "diseaseId")
            .write.parquet(evidence_dated_indirect_file)
        )


def get_association_score_by_datasource_dated(
    evidenceLink="indirect",
    excludeDatasource=[],
    diseaseId=[],
    targetId=[],
    # shuffle=False
):
    """
    Recalculate association scores by datasource over the years. Estimated running time: 2h.

    Args:
        evidenceLink (str): 'direct' or 'indirect'
        shuffle (bool):     True or False. Shuffle evidence for statistic p-value analysis

    Returns:
        Dataframe with association scores by datasource over the years. Columns:
        - diseaseId
        - targetId
        - datasourceId
        - year
        - score
    """

    # get dated evidence
    if evidenceLink == "direct":
        data = spark.read.parquet(evidence_dated_file)
        f = association_by_datasource_dated_file
    elif evidenceLink == "indirect":
        data = spark.read.parquet(evidence_dated_indirect_file)
        f = association_by_datasource_dated_indirect_file

    # exclude data_source
    if len(excludeDatasource):
        data = data.filter(~F.col("datasourceId").isin(excludeDatasource))

    if len(diseaseId):
        data = data.filter(~F.col("diseaseId").isin(diseaseId))

    if len(targetId):
        data = data.filter(~F.col("targetId").isin(targetId))

    if os.path.exists(f):
        pass

    else:

        # shuffle evd (for statistic p-value calculations)
        # if shuffle:
        #     data = (
        #         data.select(["diseaseId", "targetId"])
        #         .orderBy(F.rand())
        #         .withColumn("idx", F.monotonically_increasing_id())
        #         .join(
        #             data.select(
        #                 [
        #                     col
        #                     for col in data.columns
        #                     if (col != "diseaseId") & (col != "targetId")
        #                 ]
        #             ).withColumn("idx", F.monotonically_increasing_id()),
        #             "idx",
        #             "inner",
        #         )
        #         .drop("idx")
        #     )

        # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
        data = data.withColumn("year", F.col("year").cast(T.IntegerType())).fillna(
            last_year + 1, subset=["year"]
        )

        # get all the combinations of datasourceId vs years in the range between the first_year and the last_year set
        sourceVSyear = (
            # unique sources
            data.select("datasourceId")
            .distinct()
            .crossJoin(
                # unique years in range
                spark.createDataFrame(
                    data=[
                        [r] for r in range(first_year, last_year + 1 + 1, 1)
                    ],  # last_year+1 as a surrogate for non-dated evidence
                    schema=["year"],
                )
            )
            .repartition(400, "year")  # repartition required after crossJoin
        )

        # get all the combinations of datasourceId vs years vs disease-target score
        data = (
            sourceVSyear.join(
                data.select("diseaseId", "targetId", "datasourceId").distinct(),
                "datasourceId",
                "left",
            )
            # disease - target - datasource - year - score
            .join(
                data,
                ["diseaseId", "targetId", "datasourceId", "year"],
                "left",
            )
        )

        # prepare partition: all evidence accumulated for each disease-target-datasource-year triplet until the given year
        partition1 = (
            Window.partitionBy("diseaseId", "targetId", "datasourceId")
            .orderBy("year")
            .rangeBetween(Window.unboundedPreceding, 0)
        )

        # recalculate harmonic sum scores over the years considering the evidence accumulated
        data = (
            data.groupBy("diseaseId", "targetId", "datasourceId", "year")
            # collect list of scores FOR each year
            .agg(F.collect_list("score").alias("cum_scores"))
            # collect scores UNTIL each year
            .withColumn(
                "cum_scores", F.flatten(F.collect_list("cum_scores").over(partition1))
            )
            # remove NaNs from the cumulative scores array
            .withColumn(
                "scores_no_nan", F.expr("filter(cum_scores, x -> NOT isnan(x))")
            )
            # sort descending and take top 50 scores
            .withColumn("scores_sorted", F.reverse(F.array_sort("scores_no_nan")))
            .withColumn("top50_scores", F.expr("slice(scores_sorted, 1, 50)"))
            # generate indices (1-based) for the top 50 scores
            .withColumn("idx", F.sequence(F.lit(1), F.size("top50_scores")))
            # zip scores and indices, and divide each score by idx^2
            .withColumn(
                "weighted_scores",
                F.expr(
                    "transform(arrays_zip(top50_scores, idx), x -> x.top50_scores / pow(x.idx, 2))"
                ),
            )
            # sum the weighted scores
            .withColumn(
                "harmonic_sum",
                F.expr("aggregate(weighted_scores, 0D, (acc, x) -> acc + x)"),
            )
            # normalize by max_score
            .withColumn("harmonic_score", F.col("harmonic_sum") / F.lit(max_score))
            # select
            .select(
                "targetId",
                "diseaseId",
                "year",
                "datasourceId",
                F.col("harmonic_score").alias("score"),
            )
            # recover non-dated evidence
            .withColumn(
                "year",
                F.when(F.col("year") == last_year + 1, None).otherwise(F.col("year")),
            )
            .withColumn("sourceId", F.col("datasourceId"))
            .write.partitionBy("sourceId")
            .parquet(f)
        )


def get_association_novelty_by_datasource_dated(
    evidenceLink="indirect",
    scale=novelty_scale,
    shift=novelty_shift,
    window=novelty_window,
    diseaseId=[],
    targetId=[],
    excludeDatasource=[],
):
    """
    Calculate novelty of association scores by datasource over the years.

    Args:
        evidenceLink (str): 'direct' or 'indirect'
        scale (float):      logistic growth rate or steepness of the novelty curve
        shift (float):      x-axis value of the sigmoid's point of the novelty curve
        window (float):     range of years after the peak year to apply decay function to

    Returns:
        Dataframe with novelty of association scores by datasource over the years. Columns:
        - diseaseId
        - targetId
        - datasourceId
        - year
        - score
        - novelty
    """

    # get dated score
    if evidenceLink == "direct":
        f = association_by_datasource_dated_file
    elif evidenceLink == "indirect":
        f = association_by_datasource_dated_indirect_file
    data = spark.read.parquet(f)

    # exclude data_source
    if len(excludeDatasource):
        data = data.filter(~F.col("datasourceId").isin(excludeDatasource))

    if len(diseaseId):
        data = data.filter(~F.col("diseaseId").isin(diseaseId))

    if len(targetId):
        data = data.filter(~F.col("targetId").isin(targetId))

    # prepare partition: disease-target-datasource triplet ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId", "datasourceId").orderBy(
        "year"
    )

    # calculate novelty for disease-target-datasource over the years
    data = (
        data
        # fill non-dated evidence with last_year+1 to avoid loosing them (we'll reset them later)
        .fillna(last_year + 1, subset=["year"])
        # fill NaN score with 0 for later novelty calculation
        .fillna(0, subset=["score"])
        # for each target-disease-datasource get peaks of score shift (when current year score minus previous year score > 0)
        .select(
            "diseaseId",
            "targetId",
            "datasourceId",
            F.col("year").alias("peakYear"),
            (F.col("score") - F.lag("score", offset=1).over(partition1)).alias("peak"),
        )
        # filter peaks
        .filter(F.col("peak") > 0)
        # for each peak year, get the range of years within a window, e.g. peakYear=1991 -> range of peakYear + window=3: 1991, 1992, 1993
        .select(
            "*",
            F.posexplode(  # Returns a new row for each element with position in the given array or map. Uses the default column name pos for position, and col for elements in the array and key and value for elements in the map unless specified otherwise.
                F.sequence(  # Generate a sequence of integers from start to stop, incrementing by step. If step is not set, incrementing by 1 if start is less than or equal to stop, otherwise -1.
                    F.col("peakYear"),
                    F.col("peakYear") + window,
                )
            ).alias(
                "year-peakYear", "year"
            ),
        )
        # for each peak, calculate the novelty value at the different years within the window (novelty = peakScore/(1+exp^(scale*(year-peakYear-shift))) -> logistic function)
        # and select max. novelty value found for each year
        .groupBy("diseaseId", "targetId", "datasourceId", "year")
        .agg(
            F.round(
                F.max(
                    F.col("peak")
                    / (1 + F.exp(scale * (F.col("year-peakYear") - shift))),
                ),
                3,
            ).alias("novelty")
        )
        # add max. novelty values to original disease-target-datasource-year dataframe
        .join(
            data,
            ["diseaseId", "targetId", "datasourceId", "year"],
            "right",
        )
        # reset non-dated evidence to null year
        .replace(last_year + 1, None, subset=["year"])
        # set novelty=0 when novelty=null
        .fillna(0, subset=["novelty"])
        .withColumn("sourceId", F.col("datasourceId"))
        .write.mode("overwrite")
        .partitionBy("sourceId")
        .parquet(f)
    )


def get_association_score_by_overall_dated(
    evidenceLink="indirect",
    excludeDatasource=[],
):
    """
    Recalculate overall association scores over the years.

    Args:
        evidenceLink (str):         'direct' or 'indirect'
        excludeDatasource (list):   data_source to exclude, e.g. "chembl"

    Returns:
        Dataframe with overall association scores over the years. Columns:
        - diseaseId
        - targetId
        - year
        - score
    """

    # get dated score
    if evidenceLink == "direct":
        data = spark.read.parquet(association_by_datasource_dated_file)
        f = association_by_overall_dated_file
    elif evidenceLink == "indirect":
        data = spark.read.parquet(association_by_datasource_dated_indirect_file)
        f = association_by_overall_dated_indirect_file

    if os.path.exists(f):
        pass

    else:

        # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
        data = data.fillna(last_year + 1, subset=["year"])

        # get data_source' weights in overall score formula
        weights = spark.createDataFrame(
            data=[
                [datasourceId, str(weight)]
                for datasourceId, weight in get_weight_for_datasource()
            ],
            schema=["datasourceId", "weight"],
        )

        # exclude data_source
        if len(excludeDatasource):
            data = data.filter(~F.col("datasourceId").isin(excludeDatasource))

        # recalculate harmonic sum scores over the years considering the evidence accumulated
        data = (
            # add data_source' weights
            data.join(weights, "datasourceId", "left")
            # weight source-specific scores
            .withColumn("score", F.col("score") * F.col("weight"))
            # list of weighted scores
            .groupBy("diseaseId", "targetId", "year")
            # collect list of scores FOR each year
            .agg(F.collect_list("score").alias("cum_scores"))
            # remove NaNs from the cumulative scores array
            .withColumn(
                "scores_no_nan", F.expr("filter(cum_scores, x -> NOT isnan(x))")
            )
            # sort descending and take top 50 scores
            .withColumn("scores_sorted", F.reverse(F.array_sort("scores_no_nan")))
            .withColumn("top50_scores", F.expr("slice(scores_sorted, 1, 50)"))
            # generate indices (1-based) for the top 50 scores
            .withColumn("idx", F.sequence(F.lit(1), F.size("top50_scores")))
            # zip scores and indices, and divide each score by idx^2
            .withColumn(
                "weighted_scores",
                F.expr(
                    "transform(arrays_zip(top50_scores, idx), x -> x.top50_scores / pow(x.idx, 2))"
                ),
            )
            # sum the weighted scores
            .withColumn(
                "harmonic_sum",
                F.expr("aggregate(weighted_scores, 0D, (acc, x) -> acc + x)"),
            )
            # normalize by max_score
            .withColumn("harmonic_score", F.col("harmonic_sum") / F.lit(max_score))
            # select
            .select(
                "targetId",
                "diseaseId",
                "year",
                F.col("harmonic_score").alias("score"),
            )
            # recover non-dated evidence
            .withColumn(
                "year",
                F.when(F.col("year") == last_year + 1, None).otherwise(F.col("year")),
            )
            .write.parquet(f)
        )


def get_association_novelty_by_overall_dated(
    evidenceLink="indirect",
    scale=novelty_scale,
    shift=novelty_shift,
    window=novelty_window,
    excludeDatasource=[],
):
    """
    Calculate overall novelty of association scores over the years.

    Args:
        evidenceLink (str):         'direct'/'indirect'
        scale (float):              logistic growth rate or steepness of the novelty curve
        shift (float):              x value of the sigmoid's point of the novelty curve
        window (float):             range of years after the peak year to apply decay function to
        excludeDatasource (list):   data_source to exclude, e.g. "chembl"

    Returns:
        Dataframe with overall novelty of association scores over the years. Columns:
        - diseaseId
        - targetId
        - year
        - score
        - novelty
    """

    # get dated score
    if evidenceLink == "direct":
        f = association_by_overall_dated_file
    elif evidenceLink == "indirect":
        f = association_by_overall_dated_indirect_file
    data = spark.read.parquet(
        f
    ).persist()  # persist since we will call this dataframe more than once

    # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
    data = data.fillna(last_year + 1, subset=["year"])

    print(data.show())

    # calculate novelty for disease-target over the years
    (
        data
        # fill NaN score with 0 for later novelty calculation
        .fillna(0, subset=["score"])
        # for each target-disease get peaks of score shift (when current year score minus previous year score > 0)
        .select(
            "diseaseId",
            "targetId",
            F.col("year").alias("peakYear"),
            #  partition: disease-target pair ordered by increasing year
            (
                F.col("score")
                - F.lag("score", offset=1).over(
                    Window.partitionBy("diseaseId", "targetId").orderBy("year")
                )
            ).alias("peak"),
        )
        # filter peaks
        .filter(F.col("peak") > 0)
        # for each peak year, get the range of years within a window, e.g. peakYear=1991 -> range of peakYear + window=3: 1991, 1992, 1993
        .select(
            "*",
            F.posexplode(
                F.sequence(
                    F.col("peakYear"),
                    F.col("peakYear") + window,
                )
            ).alias("year-peakYear", "year"),
        )
        # for each peak, calculate the novelty value at the different years within the window (novelty = peakScore/(1+exp^(scale*(year-peakYear-shift))) -> logistic function)
        # and select max. novelty value found for each year
        .groupBy("diseaseId", "targetId", "year")
        .agg(
            F.round(
                F.max(
                    F.col("peak")
                    / (1 + F.exp(scale * (F.col("year-peakYear") - shift))),
                ),
                3,
            ).alias("novelty")
        )
        # add max. novelty values to original dataframe disease-target-year
        .join(
            data,
            ["diseaseId", "targetId", "year"],
            "right",
        )
        # reset non-dated evidence to null year
        .replace(last_year + 1, None, subset=["year"])
        # set novelty=0 when novelty=null
        .fillna(0, subset=["novelty"])
        .write.mode("overwrite")
        .parquet(f)
    )

    print(data.show())
    data.unpersist()


# Run novelty assessment functions
start_time = time.perf_counter()

get_indirect_evidence()

get_association_score_by_datasource_dated(
    evidenceLink="direct",
)
get_association_novelty_by_datasource_dated(
    evidenceLink="direct",
)
get_association_score_by_datasource_dated(
    evidenceLink="indirect",
)
get_association_novelty_by_datasource_dated(
    evidenceLink="indirect",
)
get_association_score_by_overall_dated(evidenceLink="direct")
get_association_novelty_by_overall_dated(evidenceLink="direct")
get_association_score_by_overall_dated(evidenceLink="indirect")
get_association_novelty_by_overall_dated(evidenceLink="indirect")


end_time = time.perf_counter()
elapsed_seconds = end_time - start_time
elapsed_td = timedelta(seconds=elapsed_seconds)
days = elapsed_td.days
hours, remainder = divmod(elapsed_td.seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"\nElapsed time: {days:02d}-{hours:02d}:{minutes:02d}:{seconds:02d}")
