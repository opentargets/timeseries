#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "26 Jun 2025"

"""
novelty_distribution.py: Distribution of novel target-disease associations in the Open Targets Platform.
"""

# Setup Google Cloud machine and sumbit job:
# gcloud dataproc clusters create cf-timeseries --image-version 2.2 --region europe-west1 --master-machine-type n1-standard-2 --secondary-worker-type spot --worker-machine-type n1-standard-4 --worker-boot-disk-size 500 --autoscaling-policy=otg-etl --optional-components=JUPYTER --enable-component-gateway --project open-targets-eu-dev
# gcloud dataproc jobs submit pyspark novelty_distribution.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"


from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window, DataFrame
from functools import reduce

# Paths
ot_platform_version = "25.03"

data_path = "gs://open-targets-data-releases/{}/output/".format(ot_platform_version)
disease_file = data_path + "disease"
target_file = data_path + "target"

results_path = "gs://ot-team/cfalaguera/{}/".format(ot_platform_version)
association_by_datasource_dated_file = results_path + "association_by_datasource_dated"
novelty_by_datasource_dated_file = results_path + "novelty_by_datasource_dated"

# Functions
prioritized_therapeutic_area = [
    "MONDO_0045024",  # "cancer or benign tumor",
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


def get_therapeutic_area_for_disease():
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
            for ranking, therapeuticArea in enumerate(prioritized_therapeutic_area)
        ],
        schema=["therapeuticArea", "ranking"],
    )

    # select top TA for each disease
    partition = Window.partitionBy("diseaseId").orderBy(F.col("ranking"))
    diseases = spark.read.parquet(disease_file).persist()
    diseases = (
        # get diseases mapped to multiple TAs
        diseases.select(
            F.col("id").alias("diseaseId"),
            F.explode("therapeuticAreas").alias("therapeuticArea"),
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


# Novels distribution
spark = SparkSession.builder.getOrCreate()

for cutoff in [0.0, 0.1, 0.5]:
    data = (
        spark.read.parquet(novelty_by_datasource_dated_file)
        .filter(
            (F.col("year") >= 2000)
            & (F.col("year") <= 2025)
            & (F.col("novelty") >= cutoff)
        )
        .join(
            get_therapeutic_area_for_disease(),
            "diseaseId",
            "inner",
        )
        # .filter(
        #     ~F.col("therapeuticArea").isin(
        #         [
        #             "GO_0008150",  # biological process
        #             "EFO_0001444",  # measurement
        #             "EFO_0002571",  # medical procedure
        #             "EFO_0000651",  # phenotype
        #             "EFO_0005932",  # animal disease
        #         ]
        #     )
        # )
        # .filter(
        #     ~F.col("diseaseId").isin(
        #         [
        #             "GO_0008150",  # biological process
        #             "EFO_0001444",  # measurement
        #             "EFO_0002571",  # medical procedure
        #             "EFO_0000651",  # phenotype
        #             "EFO_0005932",  # animal disease
        #         ]
        #     )
        # )
        .join(
            spark.read.parquet(target_file)
            # .filter(F.col("biotype") == "protein_coding")
            .select(F.col("id").alias("targetId")),
            "targetId",
            "inner",
        )
        .withColumn(
            "therapeuticAreaName",
            F.when(
                (F.col("therapeuticAreaName") == "cancer or benign tumor"),
                F.lit("Oncological"),
            )
            .when(
                (
                    F.col("therapeuticAreaName")
                    == "genetic, familial or congenital disease"
                ),
                F.lit("Congenital"),
            )
            .otherwise(F.lit("Other")),
        )
        .withColumn(
            "datasourceId",
            F.when(
                F.col("datasourceId").isin(
                    [
                        "sysbio",
                        "crispr",
                        "crispr_screen",
                        "slapenrich",
                        "gene_burden",
                        "cancer_biomarkers",
                        "reactome",
                    ]
                ),
                F.lit("other"),
            ).otherwise(F.col("datasourceId")),
        )
        .withColumn(
            "max_novelty",
            F.max("novelty").over(
                Window.partitionBy("targetId", "diseaseId", "datasourceId")
            ),
        )
        .filter(F.col("novelty") == F.col("max_novelty"))
    )

    association = data.groupby("year", "datasourceId", "therapeuticAreaName").agg(
        F.size(F.collect_set(F.concat("diseaseId", "targetId"))).alias("association")
    )

    target = (
        data.groupby("datasourceId", "targetId")
        .agg(F.min("year").alias("year"))
        .groupby("year", "datasourceId")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
    )

    # fill empty years
    (
        spark.createDataFrame(
            data=[[r] for r in range(2000, 2026, 1)],
            # data=[[r] for r in range(1900, 2026, 1)],
            schema=["year"],
        )
        .crossJoin(data.select("datasourceId").distinct())
        .crossJoin(data.select("therapeuticAreaName").distinct())
        .join(
            association,
            on=["datasourceId", "year", "therapeuticAreaName"],
            how="left",
        )
        .join(
            target,
            on=["datasourceId", "year"],
            how="left",
        )
        .fillna(0, subset=["association", "target"])
        .write.parquet(results_path + "novelty_distribution/cutoff={}".format(cutoff))
    )

    print(
        cutoff,
        data.filter(F.col("datasourceId") != "chembl")
        .select("targetId")
        .distinct()
        .count(),
    )
