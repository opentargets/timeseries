#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "10 Jul 2023"

"""
timeseries.py: Assess the evolution over time of evidence supporting disease-target associations in the Open Targets Platform. Estimated running time: 4h.
"""

# generate gcloud machine
"""
gcloud dataproc clusters create cf-timeseries --region europe-west1 --zone europe-west1-d --single-node --master-machine-type n2-highmem-128 --master-boot-disk-size 500 --image-version 2.0-debian10 --project open-targets-eu-dev
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark timeseries.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"
"""

import datetime
import os

import numpy as np
from pyspark.ml import functions as fml
from pyspark.ml.linalg import DenseVector
from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window
from pyspark.sql import types as T


# Required data
firstYear = 1970
lastYear = 2023  # datetime.date.today().year
noveltyScale = 2  # 2 for long tailed slow decays
noveltyShift = 2  # 3 for long tailed slow decays
noveltyWindow = 10
maxScore = np.sum(1 / np.arange(1, 100000 + 1) ** 2)


dataSources = [
    {
        "id": "ot_genetics_portal",
        "sectionId": "otGenetics",
        "label": "OT Genetics",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,  # needs to be a float
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#open-targets-genetics",
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
ot_platform_version = "23.06"

## data path
data_path = "gs://open-targets-data-releases/{}/output/etl/parquet/".format(
    ot_platform_version
)
evidence_file = data_path + "evidence"
diseases_file = data_path + "diseases"
literatureIndex_file = data_path + "literature/literatureIndex"

## results path
results_path = "gs://ot-team/cfalaguera/novelty/{}/".format(ot_platform_version)

evidenceIndirect_file = results_path + "evidenceIndirect"

evidenceIndirectDated_file = results_path + "evidenceIndirectDated"

evidenceDirectDated_file = results_path + "evidenceDirectDated"

associationByDatasourceIndirectOverYears_file = (
    results_path + "associationByDatasourceIndirectOverYears"
)

associationByOverallIndirectOverYears_file = (
    results_path + "associationByOverallIndirectOverYears"
)

associationByDatasourceDirectOverYears_file = (
    results_path + "associationByDatasourceDirectOverYears"
)

associationByOverallDirectOverYears_file = (
    results_path + "associationByOverallDirectOverYears"
)

associationByGeneticIndirectOverYears_file = (
    results_path + "associationByGeneticIndirectOverYears"
)
associationByLiteratureIndirectOverYears_file = (
    results_path + "associationByLiteratureIndirectOverYears"
)
associationByClinicalIndirectOverYears_file = (
    results_path + "associationByClinicalIndirectOverYears"
)

reportEvidenceIndirectDated_file = results_path + "reportEvidenceIndirectDated"
reportEvidenceDirectDated_file = results_path + "reportEvidenceDirectDated"


# Establish spark connection
spark = SparkSession.builder.getOrCreate()


# Prerequired functions
def getHarmonicScore(cumScores: DenseVector) -> float:
    """
    Calculate harmonic sum score of list of scores.

    Args:
        cumScores (DenseVector):    list of evidence score

    Returns:
        harmonicScore (float):      harmonic sum score of list of evidence score
    """

    # convert scores into numpy array (to be able to apply isnan function)
    cumScores = np.array(cumScores)

    # remove nan scores
    cumScores = cumScores[~np.isnan(cumScores)]

    # sort scores in descending order and select top50 to speed up
    cumScores = np.sort(cumScores)[::-1][:50]

    # calculate harmonic sum score = ( score1/1^2 + score2/2^2 + score3/3^2 ) / ( 1/1^2 + 1/2^2 + 1/3^2 ... up to 100000)
    harmonicScore = np.round(
        np.sum(cumScores / np.arange(1, len(cumScores) + 1) ** 2) / maxScore,
        3,
    )

    return float(harmonicScore)


def getDatasourceToWeight():
    """
    Returns list of data sources' weights for overall score.
    """

    weights = [[datasource["id"], datasource["weight"]] for datasource in dataSources]

    return weights


# Run novelty assessment functions
def getEvidence(
    evidenceLink="indirect",
    diseases_file=diseases_file,
    evidence_file=evidence_file,
):
    """
    Get evidence from OT.

    Args:
        evidenceLink (str):  'direct' or 'indirect'
        diseases_file (str): path to OT diseases file
        evidence_file (str): path to OT evidence file

    Returns:
        Dataframe with evidence. Columns:
        - targetId
        - datasourceId
        - score
        - literature
        - studyStartDate
        - diseaseId
        - drugId
        - clinicalPhase
    """

    if evidenceLink == "direct":
        return spark.read.parquet(evidence_file)

    elif evidenceLink == "indirect":
        # disease ontology expansion
        diseases = spark.read.parquet(diseases_file).select(
            F.col("id").alias("diseaseId"),
            F.explode(
                F.array_union(F.array(F.col("id")), F.col("ancestors"))
            ).alias(  # add descendant to ancestors list to keep it for later
                "specificDiseaseId"
            ),
        )

        # get indirect evidence
        evidenceIndirect = (
            spark.read.parquet(evidence_file)
            .select(
                "diseaseId",
                "targetId",
                "datasourceId",
                "score",
                "literature",
                "studyStartDate",
                "drugId",
                "clinicalPhase",
            )
            .join(diseases, "diseaseId", "inner")  # add ancestors
            .drop("diseaseId")  # drop descendants
            .withColumnRenamed("specificDiseaseId", "diseaseId")
        )

        return evidenceIndirect


def getEvidenceDated(evidenceLink="indirect"):
    """
    Map evidence to their publication year (coming from literatureIndex OT file or from evidence studyStartDate field). Estimated running time: 2 mins.

    Args:
        evidenceLink (str):             'direct' or 'indirect'

    Returns:
        Dataframe with evidence mapped to their publication year. Columns:
        - targetId
        - datasourceId
        - score
        - year
        - diseaseId
        - clinicalPhase
        - drugId
    """

    if evidenceLink == "direct":
        if os.path.exists(evidenceDirectDated_file):
            return spark.read.parquet(evidenceDirectDated_file)

    elif evidenceLink == "indirect":
        if os.path.exists(evidenceIndirectDated_file):
            return spark.read.parquet(evidenceIndirectDated_file)

    # get evidence
    evidenceDated = getEvidence(
        evidenceLink=evidenceLink
    ).persist()  #  persist since we will call this dataframe more than once

    # map evidence coming from literature to their publication year
    publications = (
        evidenceDated.filter((F.col("datasourceId") != "chembl"))
        .select(
            "diseaseId",
            "targetId",
            "datasourceId",
            "score",
            "clinicalPhase",
            "drugId",
            F.monotonically_increasing_id().alias(
                "evidenceId"
            ),  # add evidenceId to pin the original evidence for later
            F.explode_outer("literature").alias("pmid"),  # keep evidence with no pmid
        )
        # add publication year for pmids
        .join(
            spark.read.parquet(literatureIndex_file).select("pmid", "year"),
            "pmid",
            "left",  # keep evidence with no year
        )
        # get the earliest year for each evidence if there's more than one; or NULL if there's no year
        .groupBy(
            "diseaseId",
            "targetId",
            "datasourceId",
            "score",
            "clinicalPhase",
            "drugId",
            "evidenceId",
        )
        .agg(F.min("year").alias("year"))
        .drop("evidenceId")
    )

    # map evidence coming from studies to their study start year (chembl evidence)
    studies = evidenceDated.filter((F.col("datasourceId") == "chembl")).select(
        "diseaseId",
        "targetId",
        "datasourceId",
        "score",
        "clinicalPhase",
        "drugId",
        # get year from study start date
        F.split(F.col("studyStartDate"), "-").getItem(0).alias("year"),
    )

    # studies + publications
    evidenceDated = studies.unionByName(publications)

    # unpersist to free up memory
    evidenceDated.unpersist()

    return evidenceDated


def reportEvidenceDated(evidenceLink="indirect"):
    """
    Report the amount of evidence by source with or without publication date annotated. Estimated running time: 2h.

    Args:
        evidenceLink (str): 'direct' or 'indirect'

    Returns:
        Dataframe with evidence mapped to their publication year. Columns:
        - datasourceId
        - nAll
        - nDated
    """

    evidence = (
        getEvidenceDated(evidenceLink=evidenceLink)
        .groupBy("datasourceId")
        .agg(F.count("*").alias("nAll"))
        .join(
            getEvidenceDated(evidenceLink=evidenceLink)
            .filter(F.col("year").isNotNull())
            .groupBy("datasourceId")
            .agg(F.count("*").alias("nDated")),
            "datasourceId",
            "left",
        )
    )

    return evidence


def getScoreByDatasourceOverYears(datasourceId, evidenceLink="indirect"):
    """
    Recalculate association scores by datasource over the years.

    Args:
        evidenceLink (str): 'direct' or 'indirect'

    Returns:
        Dataframe with association scores by datasource over the years. Columns:
        - diseaseId
        - targetId
        - datasourceId
        - year
        - score
    """

    # get dated evidence
    scoreByDatasourceOverYears = getEvidenceDated(evidenceLink=evidenceLink).filter(
        F.col("datasourceId") == datasourceId
    )

    # fill non-dated evidence with lastYear + 1 to avoid loosing them (we'll reset them later)
    scoreByDatasourceOverYears = scoreByDatasourceOverYears.fillna(
        str(lastYear + 1), subset=["year"]  # it needs to be a str to avoid errors
    )

    # get all the combinations of datasourceId vs years in the range between the firstYear and the lastYear set
    sourceVSyear = (
        # unique sources
        scoreByDatasourceOverYears.select("datasourceId")
        .distinct()
        .crossJoin(
            # unique years in range
            spark.createDataFrame(
                data=[
                    [r] for r in range(firstYear, lastYear + 1 + 1, 1)
                ],  # lastYear + 1 as a surrogate for non-dated evidence
                schema=["year"],
            )
        )
        .repartition(400, "year")  # repartition required after crossJoin
    )

    # get all the combinations of datasourceId vs years vs disease-target score
    scoreByDatasourceOverYears = (
        # disease - target - datasource
        scoreByDatasourceOverYears.select(
            "diseaseId", "targetId", "datasourceId"
        ).distinct()
        # disease - target - datasource - year
        .join(
            F.broadcast(sourceVSyear), "datasourceId", "right"
        )  # broadcast is applied by default to <100 MB dataframes, in this case, force broadcast
        # disease - target - datasource - year - score
        .join(
            scoreByDatasourceOverYears,
            ["diseaseId", "targetId", "datasourceId", "year"],
            "left",
        )
    )

    # recalculate harmonic sum scores over the years considering the evidence accumulated

    # register udf function
    harmonicScore = F.udf(getHarmonicScore, T.DoubleType())

    # prepare partition: all evidence accumulated for each disease-target-datasource triplet until the given year
    partition1 = (
        Window.partitionBy("diseaseId", "targetId", "datasourceId")
        .orderBy("year")
        .rangeBetween(
            Window.unboundedPreceding, 0
        )  # 0 = current year; unboundedPreceding = years prior to current year
    )

    # apply functions to partition: for each target-disease-datasourece-year calculate harmonic sum score
    scoreByDatasourceOverYears = (
        scoreByDatasourceOverYears.select(
            "diseaseId",
            "targetId",
            "datasourceId",
            "year",
            harmonicScore(
                fml.array_to_vector(F.collect_list("score").over(partition1))
            ).alias("score"),
        )
        .distinct()  # get rid of rows multiplicity coming from different original scores
        .replace(
            float("nan"), None, subset=["score"]
        )  # convert score = nan values into null values for later (pyspark filter function misshandles nan values)
        # reset non-dated evidence to null year
        .replace(
            str(lastYear + 1), None, subset=["year"]
        )  # it needs to be a string to avoid errors
    )

    return scoreByDatasourceOverYears


def getNoveltyByDatasourceOverYears(
    datasourceId,
    evidenceLink="indirect",
    scale=noveltyScale,
    shift=noveltyShift,
    window=noveltyWindow,
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

    scoreByDatasourceOverYears = getScoreByDatasourceOverYears(
        evidenceLink=evidenceLink, datasourceId=datasourceId
    ).persist()  #  persist since we will call this dataframe more than once

    # prepare partition: disease-target-datasource triplet ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId", "datasourceId").orderBy(
        "year"
    )

    # calculate novelty for disease-target-datasource over the years
    noveltyByDatasourceOverYears = (
        scoreByDatasourceOverYears
        # fill non-dated evidence with lastYear + 1 to avoid loosing them (we'll reset them later)
        .fillna(lastYear + 1, subset=["year"])
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
            F.posexplode(
                F.sequence(
                    F.col("peakYear"),
                    F.col("peakYear") + window,
                )
            ).alias("year-peakYear", "year"),
        )
        # for each peak, calculate the novelty value at the different years within the window (novelty = peakScore/(1+exp^(scale*(year-peakYear-shift))) -> logistic function)
        # and select max. novelty value found for each year
        .groupBy("diseaseId", "targetId", "datasourceId", "year").agg(
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
            scoreByDatasourceOverYears,
            ["diseaseId", "targetId", "datasourceId", "year"],
            "right",
        )
        # reset non-dated evidence to null year
        .replace(lastYear + 1, None, subset=["year"])
        # set novelty=0 when novelty=null
        .fillna(0, subset=["novelty"])
    )

    # scoreByDatasourceOverYears.unpersist()

    return noveltyByDatasourceOverYears


def getScoreByOverallOverYears(
    evidenceLink="indirect",
    excludeDatasource=[],
):
    """
    Recalculate overall association scores over the years.

    Args:
        evidenceLink (str):         'direct' or 'indirect'
        excludeDatasource (list):   datasources to exclude, e.g. "chembl"

    Returns:
        Dataframe with overall association scores over the years. Columns:
        - diseaseId
        - targetId
        - year
        - score
    """

    # get association score by datasource over the years
    if evidenceLink == "direct":
        scoreByDatasourceOverYears = spark.read.parquet(
            associationByDatasourceDirectOverYears_file
        ).persist()  #  persist since we will call this dataframe more than once
    elif evidenceLink == "indirect":
        scoreByDatasourceOverYears = spark.read.parquet(
            associationByDatasourceIndirectOverYears_file
        ).persist()  #  persist since we will call this dataframe more than once

    # fill non-dated evidence with lastYear + 1 to avoid loosing them (we'll reset them later)
    scoreByDatasourceOverYears = scoreByDatasourceOverYears.fillna(
        lastYear + 1, ["year"]
    )

    # get datasources' weights in overall score formula
    weights = spark.createDataFrame(
        data=[
            [datasourceId, str(weight)]
            for datasourceId, weight in getDatasourceToWeight()
        ],
        schema=["datasourceId", "weight"],
    )

    # exclude datasources
    if len(excludeDatasource):
        scoreByDatasourceOverYears = scoreByDatasourceOverYears.filter(
            ~F.col("datasourceId").isin(excludeDatasource)
        )

    # recalculate harmonic sum scores over the years considering the evidence accumulated

    # register udfs
    harmonicScore = F.udf(getHarmonicScore, T.DoubleType())

    # prepare partition: all evidence accumulated for each disease-target pair until the given year
    partition1 = Window.partitionBy("diseaseId", "targetId", "year")

    # apply functions to partitions: for each target-disease-year calculate overall harmonic sum score
    scoreByOverallOverYears = (
        # add datasources' weights
        scoreByDatasourceOverYears.join(weights, "datasourceId", "left")
        # calculate overall harmonic score
        .select(
            "diseaseId",
            "targetId",
            "year",
            harmonicScore(
                fml.array_to_vector(  # vector makes the following calculations faster
                    F.collect_list(F.col("score") * F.col("weight")).over(partition1)
                )
            ).alias("score"),
            # get rid of "score" and "weight" multiplicity
        ).distinct()
        # convert score = nan values into null values (pyspark filter function misshandles nan values)
        .replace(float("nan"), None, subset=["score"])
        # reset non-dated evidence to null year
        .replace(lastYear + 1, None, subset=["year"])
    )

    scoreByDatasourceOverYears.unpersist()

    return scoreByOverallOverYears


def getNoveltyByOverallOverYears(
    evidenceLink="indirect",
    scale=noveltyScale,
    shift=noveltyShift,
    window=noveltyWindow,
    excludeDatasource=[],
):
    """
    Calculate overall novelty of association scores over the years.

    Args:
        evidenceLink (str):         'direct'/'indirect'
        scale (float):              logistic growth rate or steepness of the novelty curve
        shift (float):              x value of the sigmoid's point of the novelty curve
        window (float):             range of years after the peak year to apply decay function to
        excludeDatasource (list):   datasources to exclude, e.g. "chembl"

    Returns:
        Dataframe with overall novelty of association scores over the years. Columns:
        - diseaseId
        - targetId
        - year
        - score
        - novelty
    """

    scoreByOverallOverYears = (
        getScoreByOverallOverYears(
            evidenceLink=evidenceLink, excludeDatasource=excludeDatasource
        )
        # fill non-dated evidence with lastYear + 1 to avoid loosing them (we'll reset them later)
        .fillna(lastYear + 1, ["year"])
        # persist since we will call this dataframe more than once
        .persist()
    )

    # prepare partition: disease-target pair ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId").orderBy("year")

    # calculate novelty for disease-target over the years
    noveltyByOverallOverYears = (
        scoreByOverallOverYears
        # fill NaN score with 0 for later novelty calculation
        .fillna(0, subset=["score"])
        # for each target-disease get peaks of score shift (when current year score minus previous year score > 0)
        .select(
            "diseaseId",
            "targetId",
            F.col("year").alias("peakYear"),
            (F.col("score") - F.lag("score", offset=1).over(partition1)).alias("peak"),
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
        .groupBy("diseaseId", "targetId", "year").agg(
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
            scoreByOverallOverYears,
            ["diseaseId", "targetId", "year"],
            "right",
        )
        # reset non-dated evidence to null year
        .replace(lastYear + 1, None, subset=["year"])
        # set novelty=0 when novelty=null
        .fillna(0, subset=["novelty"])
    )

    scoreByOverallOverYears.unpersist()

    return noveltyByOverallOverYears


def writeParquet(dataframe, filename):
    """
    Write dataframe into parquet file.

    Args:
        dataframe (dataframe):  dataframe
        filename (str):         output file name

    Returns:
        Parquet file with dataframe

    """

    if not os.path.exists(filename):
        print("writting {}...".format(filename))
        dataframe.write.parquet(filename)
        print("{} succesfully generated! :)".format(filename))
        spark.catalog.clearCache()  # remove all tables to free up space


if 0:
    # dating evidence
    writeParquet(
        dataframe=getEvidenceDated(evidenceLink="indirect"),
        filename=evidenceIndirectDated_file,
    )

    writeParquet(
        dataframe=reportEvidenceDated(evidenceLink="indirect"),
        filename=reportEvidenceIndirectDated_file,
    )

    # datasource novelty
    for datasourceId in [datasource["id"] for datasource in dataSources]:
        writeParquet(
            dataframe=getNoveltyByDatasourceOverYears(
                evidenceLink="indirect", datasourceId=datasourceId
            ),
            filename=associationByDatasourceIndirectOverYears_file
            + "/sourceId="
            + datasourceId,
        )

    # overall novelty
    writeParquet(
        dataframe=getNoveltyByOverallOverYears(evidenceLink="indirect"),
        filename=associationByOverallIndirectOverYears_file,
    )

    # genetic novelty
    writeParquet(
        dataframe=getNoveltyByOverallOverYears(
            evidenceLink="indirect", excludeDatasource=["chembl", "europepmc"]
        ),
        filename=associationByGeneticIndirectOverYears_file,
    )

    # literature novelty
    writeParquet(
        dataframe=getNoveltyByOverallOverYears(
            evidenceLink="indirect",
            excludeDatasource=[
                datasource["id"]
                for datasource in dataSources
                if datasource["id"] != "europepmc"
            ],
        ),
        filename=associationByLiteratureIndirectOverYears_file,
    )

    # clinical novelty
    writeParquet(
        dataframe=getNoveltyByOverallOverYears(
            evidenceLink="indirect",
            excludeDatasource=[
                datasource["id"]
                for datasource in dataSources
                if datasource["id"] != "chembl"
            ],
        ),
        filename=associationByClinicalIndirectOverYears_file,
    )
