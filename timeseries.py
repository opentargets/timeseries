#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "24 May 2023"

"""
timeseries.py: Assess the evolution over time of evidence supporting disease-target associations in the Open Targets Platform. Estimated running time: 2h.
"""

# generate gcloud machine
"""
gcloud dataproc clusters create cf-novelty --region europe-west1 --zone europe-west1-d --single-node --master-machine-type n2-standard-80 --master-boot-disk-size 500 --image-version 2.0-debian10 --project open-targets-eu-devgcloud dataproc clusters create cf-novelty-highmem --region europe-southwest1 --zone europe-southwest1-b --single-node --master-machine-type n2-standard-64 --master-boot-disk-size 2000 --image-version 2.0-debian10 --project open-targets-eu-dev
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark novelty.py --cluster=cf-novelty --project=open-targets-eu-dev --region="europe-west1"
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
firstYear = 1950
lastYear = 2023  # datetime.date.today().year
noveltyScale = 2  # 2 for long tailed slow decays
noveltyShift = 2  # 3 for long tailed slow decays
noveltyWindow = 10
maxScore = np.sum(1 / np.arange(1, 100000) ** 2)
noveltyCutoff = 0.1


# Paths
ot_platform_version = "/23.02/"

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

associationByDatasourceIndirectOverYearsSignature_file = (
    results_path + "associationByDatasourceIndirectOverYearsSignature"
)


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

    weights = [
        ["cancer_biomarkers", 0.5],
        ["cancer_gene_census", 1],
        ["chembl", 1],
        ["clingen", 1],
        ["crispr", 1],
        ["encore", 0.5],
        ["europepmc", 0.2],
        ["eva", 1],  # clinvar
        ["eva_somatic", 1],  # clinvar (somatic)
        ["expression_atlas", 0.2],
        ["gene2phenotype", 1],
        ["gene_burden", 1],
        ["genomics_england", 1],
        ["impc", 0.2],
        ["intogen", 1],
        ["orphanet", 1],
        ["ot_crispr", 0.5],
        ["ot_crispr_validation", 0.5],
        ["ot_genetics_portal", 1],
        ["progeny", 0.5],
        ["reactome", 1],
        ["slapenrich", 0.5],
        ["sysbio", 0.5],
        ["uniprot_literature", 1],
        ["uniprot_variants", 1],
    ]

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
    Map evidence to their publication year (coming from literatureIndex OT file or from evidence studyStartDate field).

    Args:
        evidenceLink (str): 'direct' or 'indirect'

    Returns:
        Dataframe with evidence mapped to their publication year. Columns:
        - targetId
        - datasourceId
        - score
        - year
        - diseaseId
    """

    # get evidence
    evidenceDated = getEvidence(evidenceLink=evidenceLink).persist()

    # map evidence coming from literature to their publicaition year
    publications = (
        evidenceDated.select(
            "diseaseId",
            "targetId",
            "datasourceId",
            "score",
            F.monotonically_increasing_id().alias(
                "evidenceId"
            ),  # add evidenceId to pin the original evidence for later
            F.explode("literature").alias("pmid"),
        )
        # add publication year for pmids
        .join(
            spark.read.parquet(literatureIndex_file)
            .select("pmid", "year")
            .distinct()
            .filter(F.col("year").isNotNull()),  # filter out evidence without a year
            "pmid",
            "inner",
        )
        # get earliest year for each evidence (if there's more than one)
        .groupBy("diseaseId", "targetId", "datasourceId", "score", "evidenceId")
        .agg(F.min("year").alias("year"))
        .drop("evidenceId")
    )

    # map evidence coming from studies to their study start year (chembl evidence)
    studies = evidenceDated.filter(
        (F.col("datasourceId") == "chembl")
        & (F.col("studyStartDate").isNotNull())  # filter out evidence without a year
    ).select(
        "diseaseId",
        "targetId",
        "datasourceId",
        "score",
        # get year from study start date
        F.split(F.col("studyStartDate"), "-").getItem(0).alias("year"),
    )

    # studies + publications
    evidenceDated = studies.unionByName(publications)

    return evidenceDated


def getScoreByDatasourceOverYears(evidenceLink="indirect"):
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

    # get evidence with year
    scoreByDatasourceOverYears = getEvidenceDated(evidenceLink=evidenceLink).persist()

    # get all the combinations of datasourceId vs years in the range between the firstYear and the lastYear set
    sourceVSyear = (
        # unique sources
        scoreByDatasourceOverYears.select("datasourceId")
        .distinct()
        .crossJoin(
            # unique years in range
            spark.createDataFrame(
                data=[[r] for r in range(firstYear, lastYear + 1, 1)],
                schema=["year"],
            )
        )
        .repartition(400, "datasourceId")  # repartition required after crossJoin
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
            float("nan"), None
        )  # convert score = nan values into null values for later (pyspark filter function misshandles nan values)
    )

    return scoreByDatasourceOverYears


def getNoveltyByDatasourceOverYears(
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

    # get association score by datasource over the years
    if evidenceLink == "direct":
        if os.path.exists(associationByDatasourceDirectOverYears_file):
            return spark.read.parquet(associationByDatasourceDirectOverYears_file)
        else:
            noveltyByDatasourceOverYears = getScoreByDatasourceOverYears(
                evidenceLink=evidenceLink
            ).persist()
    elif evidenceLink == "indirect":
        if os.path.exists(associationByDatasourceIndirectOverYears_file):
            return spark.read.parquet(associationByDatasourceIndirectOverYears_file)
        else:
            noveltyByDatasourceOverYears = getScoreByDatasourceOverYears(
                evidenceLink=evidenceLink
            ).persist()

    # prepare partition: disease-target-datasource triplet ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId", "datasourceId").orderBy(
        "year"
    )

    # calculate novelty for disease-target-datasource over the years
    noveltyByDatasourceOverYears = (
        noveltyByDatasourceOverYears.fillna(
            0, subset=["score"]
        )  # fill NaN score with 0 for later novelty calculation
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
            noveltyByDatasourceOverYears,
            ["diseaseId", "targetId", "datasourceId", "year"],
            "right",
        )
    )

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
        if os.path.exists(associationByDatasourceDirectOverYears_file):
            scoreByOverallOverYears = spark.read.parquet(
                associationByDatasourceDirectOverYears_file
            ).persist()  #  persist since we will call this dataframe more than once
        else:
            scoreByOverallOverYears = getScoreByDatasourceOverYears(
                evidenceLink=evidenceLink
            ).persist()

    elif evidenceLink == "indirect":
        if os.path.exists(associationByDatasourceIndirectOverYears_file):
            scoreByOverallOverYears = spark.read.parquet(
                associationByDatasourceIndirectOverYears_file
            ).persist()  #  persist since we will call this dataframe more than once
        else:
            scoreByOverallOverYears = getScoreByDatasourceOverYears(
                evidenceLink=evidenceLink
            ).persist()

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
        scoreByOverallOverYears = scoreByOverallOverYears.filter(
            ~F.col("datasourceId").isin(excludeDatasource)
        )

    # recalculate harmonic sum scores over the years considering the evidence accumulated

    # register udfs
    harmonicScore = F.udf(getHarmonicScore, T.DoubleType())

    # prepare partition: all evidence accumulated for each disease-target pair until the given year
    partition1 = Window.partitionBy("diseaseId", "targetId", "year")

    # apply functions to partitions: for each target-disease-year calculate overall harmonic sum score
    scoreByOverallOverYears = (
        scoreByOverallOverYears.join(
            weights, "datasourceId", "left"
        )  # add datasources' weights
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
        )
        .distinct()  # get rid of "score" and "weight" multiplicity
        .replace(
            float("nan"), None
        )  # convert score = nan values into null values (pyspark filter function misshandles nan values)
    )

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

    # get association score by datasource over the years
    if evidenceLink == "direct":
        if os.path.exists(associationByOverallDirectOverYears_file):
            return spark.read.parquet(associationByOverallDirectOverYears_file)
        else:
            noveltyByOverallOverYears = getScoreByOverallOverYears(
                evidenceLink=evidenceLink, excludeDatasource=excludeDatasource
            ).persist()

    elif evidenceLink == "indirect":
        if os.path.exists(associationByOverallIndirectOverYears_file):
            return spark.read.parquet(associationByOverallIndirectOverYears_file)
        else:
            noveltyByOverallOverYears = getScoreByOverallOverYears(
                evidenceLink=evidenceLink, excludeDatasource=excludeDatasource
            ).persist()

    # prepare partition: disease-target pair ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId").orderBy("year")

    # calculate novelty for disease-target over the years
    noveltyByOverallOverYears = (
        noveltyByOverallOverYears.fillna(
            0, subset=["score"]
        )  # fill NaN score with 0 for later novelty calculation
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
            noveltyByOverallOverYears,
            ["diseaseId", "targetId", "year"],
            "right",
        )
    )

    return noveltyByOverallOverYears


def getNoveltySignatureByDatasource(
    evidenceLink="indirect", noveltyCutoff=noveltyCutoff
):
    """
    Generate novelty signatures representing associations' novelty evolution.
    Novelty signatures consist of a binary vector with values equal 1 for those years when the novelty surpasess the novelty cutoff.

    Args:
        evidenceLink (str):         'direct'/'indirect'
        noveltyCutoff (float):      novelty value cutoff required to set a 1 in the binary vector.

    Returns:
        Dataframe with by-datasource novelty signatures of associations. Columns:
        - diseaseId
        - targetId
        - datasourceId
        - signature
    """

    # get novelty by datasource over the years
    if evidenceLink == "direct":
        signature = spark.read.parquet(associationByDatasourceDirectOverYears_file)
    elif evidenceLink == "indirect":
        signature = spark.read.parquet(associationByDatasourceIndirectOverYears_file)

    # prepare partition: disease-target-datasourc triplet ordered by increasing year
    partition1 = Window.partitionBy("diseaseId", "targetId", "datasourceId").orderBy(
        "year"
    )

    # get novelty by datasource over the years
    signature = (
        signature
        # convert novelty values into 1/0
        .withColumn(
            "binNovelty", F.when(F.col("novelty") >= noveltyCutoff, 1).otherwise(0)
        )
        # get novelty signature
        .select(
            "diseaseId",
            "targetId",
            "datasourceId",
            "year",
            F.collect_list("binNovelty").over(partition1).alias("noveltySignature"),
        )
        .filter(F.col("year") == 2023)  # remove multiplicity
        .drop("year")
        .filter(
            F.array_contains("noveltySignature", 1)
        )  # consider only those datasources with non-only-zeros signatures (novelty values above the cutoff)
    )

    return signature


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


if 0:
    writeParquet(
        dataframe=getEvidence("indirect"),
        filename=evidenceIndirect_file,
    )

    writeParquet(
        dataframe=getNoveltyByDatasourceOverYears(evidenceLink="indirect"),
        filename=associationByDatasourceIndirectOverYears_file,
    )

    writeParquet(
        dataframe=getNoveltyByOverallOverYears(evidenceLink="indirect"),
        filename=associationByOverallIndirectOverYears_file,
    )

    writeParquet(
        dataframe=getNoveltySignatureByDatasource(evidenceLink="indirect"),
        filename=associationByDatasourceIndirectOverYearsSignature_file,
    )
