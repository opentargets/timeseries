#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "6 Jul 2023"

"""
novels.py: List novel associations in the Open Targets Platform.
"""

# generate gcloud machine
"""
gcloud dataproc clusters create cf-novelty --region europe-west1 --zone europe-west1-d --single-node --master-machine-type n2-standard-80 --master-boot-disk-size 500 --image-version 2.0-debian10 --project open-targets-eu-dev
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark novels.py --cluster=cf-novelty --project=open-targets-eu-dev --region="europe-west1"
"""

import os

from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window

# Required data
year = 2023
prioritizedTherapeuticArea = [
    "MONDO_0045024",  # "cell proliferation disorder",
    "EFO_0005741",  # "infectious disease",
    "OTAR_0000014",  # "pregnancy or perinatal disease",
    "EFO_0005932",  # "animal disease",
    "MONDO_0024458",  # "disease of visual system",
    "EFO_0000319",  # "cardiovascular disease",
    "EFO_0009605",  # "pancreas disease",
    "EFO_0010282",  # "gastrointestinal disease",
    "OTAR_0000017",  # "reproductive system or breast disease",
    "EFO_0010285",  # "integumentary system disease",
    "EFO_0001379",  # "endocrine system disease",
    "OTAR_0000010",  # "respiratory or thoracic disease",
    "EFO_0009690",  # "urinary system disease",
    "OTAR_0000006",  # "musculoskeletal or connective tissue disease",
    "MONDO_0021205",  # "disease of ear",
    "EFO_0000540",  # "immune system disease",
    "EFO_0005803",  # "hematologic disease",
    "EFO_0000618",  # "nervous system disease",
    "MONDO_0002025",  # "psychiatric disorder",
    "MONDO_0024297",  # "nutritional or metabolic disease",
    "OTAR_0000018",  # "genetic, familial or congenital disease",
    "OTAR_0000009",  # "injury, poisoning or other complication",
    "EFO_0000651",  # "phenotype",
    "EFO_0001444",  # "measurement",
    "GO_0008150",  # "biological process"
]
noveltyCutoff = 0.1

# Paths
ot_platform_version = "23.06"

## data path
data_path = "gs://open-targets-data-releases/{}/output/etl/parquet/".format(
    ot_platform_version
)
diseases_file = data_path + "diseases"
targets_file = data_path + "targets"
associationByOverallDirect_file = data_path + "associationByOverallDirect"
targetsPrioritisation_file = (
    data_path + "targetsPriorisation"
)  # remember to change this to "targetPrioritisation" when the data team fixes the typo

## results path
results_path = "gs://ot-team/cfalaguera/novelty/{}/".format(ot_platform_version)

associationByDatasourceIndirectOverYears_file = (
    results_path + "associationByDatasourceIndirectOverYears"
)

associationByOverallIndirectOverYears_file = (
    results_path + "associationByOverallIndirectOverYears"
)

novels_file = results_path + "novels"

# Establish spark connection
spark = SparkSession.builder.getOrCreate()


def getTherapeuticAreaForDisease():
    """
    Map diseases to their prioritized therapeutic area.

    Returns:
        pyspark dataframe with the following columns:
        - diseaseId
        - therapeuticArea
        - therapeuticAreaName
    """

    # get TAs ranked by priority
    therapeuticAreas = spark.createDataFrame(
        data=[
            [therapeuticArea, ranking]
            for ranking, therapeuticArea in enumerate(
                # otPlatform.prioritizedTherapeuticArea
                prioritizedTherapeuticArea
            )
        ],
        schema=["therapeuticArea", "ranking"],
    )

    # select top TA for each disease
    partition = Window.partitionBy("diseaseId").orderBy(F.col("ranking"))
    diseases = spark.read.parquet(diseases_file).persist()
    diseases = (
        # get diseases mapped to multiple TAs
        diseases.select(
            F.col("id").alias("diseaseId"),
            F.explode_outer("therapeuticAreas").alias("therapeuticArea"),
        )
        # add TA ranking
        .join(therapeuticAreas, "therapeuticArea", "left")
        # fill with 1000 those TA not included in the ranking (to avoid losing them)
        .fillna(1000, subset=["ranking"])
        # select top ranking TA for each disease
        .withColumn("row", F.row_number().over(partition))
        .filter(F.col("row") == 1)
        .drop("row", "ranking")
        # add TA name
        .join(
            diseases.select(
                F.col("id").alias("therapeuticArea"),
                F.col("name").alias("therapeuticAreaName"),
            ),
            "therapeuticArea",
            "left",
        )
    )

    return diseases


def getNovels():
    """
    List novel target-disease associations in the current year and Platform version.
    These are target-disease associations that fulfil ALL the following conditions:
        i) Have at least one evidence linking them directly (to avoid related diseases redundancy)
        ii) Have an overall novelty value >= noveltyCutoff in current year
        iii) It's the first time in a window of -2 years that the novelty crosses the noveltyCutoff (both overall and by-datasource novelties)
        iv) Their score values when including or not non-dated evidence are the same
        v) Have never been (pre-)clinically explored, neither in the past nor now


    Returns:
        Parquet file with list of novels
    """

    novels = (
        # i) Have at least one evidence linking them directly (to avoid related diseases redundancy)v
        spark.read.parquet(associationByOverallDirect_file)
        .select("diseaseId", "targetId")
        # ii) Have an overall novelty value >= noveltyCutoff in current year
        .join(
            spark.read.parquet(associationByOverallIndirectOverYears_file).filter(
                (F.col("year") == year) & (F.col("novelty") >= noveltyCutoff)
            ),
            ["diseaseId", "targetId"],
            "inner",
        )
        # iii) It's the first time in a window of -2 years that the novelty crosses the noveltyCutoff (both overall and by-datasource novelties)
        .join(
            spark.read.parquet(associationByDatasourceIndirectOverYears_file)
            .filter((F.col("year") < year - 2) & (F.col("novelty") >= noveltyCutoff))
            .select("diseaseId", "targetId"),
            ["diseaseId", "targetId"],
            "anti",
        )
        # iv) Their score values when including or not non-dated evidence are the same
        .join(
            spark.read.parquet(associationByOverallIndirectOverYears_file)
            .filter((F.col("year").isNull()))
            .select("diseaseId", "targetId", F.col("score").alias("nullScore")),
            ["diseaseId", "targetId"],
            "left",
        )
        .filter(F.col("nullScore") == F.col("score"))
        # v) Have never been (pre-)clinically explored, neither in the past nor now
        .join(
            spark.read.parquet(associationByDatasourceIndirectOverYears_file)
            .filter(
                (F.col("year").isNull())
                & (F.col("datasourceId") == "chembl")
                & (F.col("score") > 0)
            )
            .select("diseaseId", "targetId"),
            ["diseaseId", "targetId"],
            "anti",
        )
        # add TA
        .join(getTherapeuticAreaForDisease(), "diseaseId", "left")
        # add supporting datasources
        .join(
            spark.read.parquet(associationByDatasourceIndirectOverYears_file)
            .filter((F.col("year") == year) & (F.col("novelty") >= noveltyCutoff))
            .select("diseaseId", "targetId", "datasourceId"),
            ["diseaseId", "targetId"],
            "left",
        )
        .groupBy(
            "therapeuticArea",
            "targetId",
            "diseaseId",
            "novelty",
            "score",
        )
        .agg(
            F.concat_ws(", ", F.sort_array(F.collect_list("datasourceId"))).alias(
                "supportingDatasources"
            )
        )
        # add disease name
        .join(
            spark.read.parquet(diseases_file).select(
                F.col("id").alias("diseaseId"), F.col("name").alias("diseaseName")
            ),
            "diseaseId",
            "left",
        )
        # add target symbol
        .join(
            spark.read.parquet(targets_file).select(
                F.col("id").alias("targetId"),
                F.col("approvedSymbol").alias("targetSymbol"),
            ),
            "targetId",
            "left",
        )
        # add therapeutic area name
        .join(
            spark.read.parquet(diseases_file).select(
                F.col("id").alias("therapeuticArea"),
                F.col("name").alias("therapeuticAreaName"),
            ),
            "therapeuticArea",
            "left",
        )
        # order
        .select(
            "therapeuticAreaName",
            "diseaseName",
            "targetSymbol",
            "novelty",
            "score",
            "supportingDatasources",
            "therapeuticArea",
            "diseaseId",
            "targetId",
        )
        # add target prioritisation
        .join(
            spark.read.parquet(targetsPrioritisation_file),
            "targetId",
            "left",
        )
        # order
        .orderBy("novelty", "targetId", ascending=False)
    )

    return novels


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


if 1:
    writeParquet(dataframe=getNovels(), filename=novels_file)
