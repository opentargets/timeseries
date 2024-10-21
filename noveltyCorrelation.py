#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "10 May 2024"

# generate gcloud machine
"""
noveltyCorrelation.py: Analyse the co-ocurrance of novelty peaks across resources of the same type.
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark noveltyCorrelation.py --cluster=cf-timeseries2 --project=open-targets-eu-dev --region="europe-west1"
"""

import sys
import itertools
from functools import reduce
import numpy as np
from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window, DataFrame
from scipy import stats
import pandas as pd

correlationNovelty_file = (
    "gs://ot-team/cfalaguera/correlationNovelty/correlationNovelty"
)

spark = SparkSession.builder.getOrCreate()

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
excludeTherapeuticArea = [
    "GO_0008150",  # biological process
    "EFO_0001444",  # measurement
    "EFO_0002571",  # medical procedure
]


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
    diseases = spark.read.parquet(
        "gs://open-targets-data-releases/23.06/output/etl/parquet/diseases"
    ).persist()
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


# pairs
if 1:
    for evidenceLink in (
        "direct",
        # "indirect"
    ):
        print(evidenceLink)
        if evidenceLink == "direct":
            associations = "gs://ot-team/cfalaguera/novelty/23.06/associationByDatasourceDirectOverYears"
        elif evidenceLink == "indirect":
            associations = "gs://ot-team/cfalaguera/novelty/23.06/associationByDatasourceIndirectOverYears"

        results2 = {}
        for dataset in ["Shuffled", ""]:
            print(dataset)
            data = (
                spark.read.parquet(associations + dataset)
                .select("diseaseId", "targetId", "datasourceId", "year", "novelty")
                .filter(F.col("year").isNotNull() & (F.col("novelty") > 0))
                # only protein_coding targets
                .join(
                    spark.read.parquet(
                        "gs://open-targets-data-releases/23.06/output/etl/parquet/targets"
                    )
                    .filter(F.col("biotype") == "protein_coding")
                    .select(F.col("id").alias("targetId")),
                    "targetId",
                    "inner",
                )
                # exclude TA
                .join(
                    getTherapeuticAreaForDisease().select(
                        "diseaseId", "therapeuticArea", "therapeuticAreaName"
                    ),
                    "diseaseId",
                    "left",
                )
                .filter(~F.col("therapeuticArea").isin(excludeTherapeuticArea))
                .filter(~F.col("diseaseId").isin(prioritizedTherapeuticArea))
                # filter clinvar in 2013 (noisy)
                # .filter(
                #     (F.col("datasourceId") != "eva")
                #     | ((F.col("datasourceId") == "eva") & (F.col("year") > 2015))
                # )
                # max. novelty by datasource
                .withColumn(
                    "maxNovelty",
                    F.max("novelty").over(
                        Window.partitionBy("diseaseId", "targetId", "datasourceId")
                    ),
                )
                .filter(F.col("novelty") == F.col("maxNovelty"))
                .groupby("diseaseId", "targetId", "datasourceId", "maxNovelty")
                .agg(F.min("year").alias("year"))
            ).persist()

            datasources = (
                data.select("datasourceId").distinct().toPandas().datasourceId.tolist()
            )

            results = []
            for datasourceIdA, datasourceIdB in itertools.product(
                datasources,
                repeat=2,
            ):
                results.append(
                    data.filter(F.col("datasourceId") == datasourceIdA)
                    .select(
                        "diseaseId",
                        "targetId",
                        F.col("datasourceId").alias("datasourceIdA"),
                        F.col("maxNovelty").alias("maxNoveltyA"),
                        F.col("year").alias("yearA"),
                    )
                    .join(
                        data.filter(F.col("datasourceId") == datasourceIdB).select(
                            "diseaseId",
                            "targetId",
                            F.col("datasourceId").alias("datasourceIdB"),
                            F.col("maxNovelty").alias("maxNoveltyB"),
                            F.col("year").alias("yearB"),
                        ),
                        ["diseaseId", "targetId"],
                        "inner",
                    )
                )

            results = reduce(DataFrame.unionByName, results).repartition(
                400, "datasourceIdA"
            )

            results2[dataset] = results

        (
            results2[""]
            .select(
                F.col("yearA").alias("realYearA"),
                F.col("yearB").alias("realYearB"),
                F.col("maxNoveltyA").alias("realMaxNoveltyA"),
                F.col("maxNoveltyB").alias("realMaxNoveltyB"),
                "diseaseId",
                "targetId",
                "datasourceIdA",
                "datasourceIdB",
            )
            .join(
                results2["Shuffled"].select(
                    F.col("yearA").alias("shuffledYearA"),
                    F.col("yearB").alias("shuffledYearB"),
                    F.col("maxNoveltyA").alias("shuffledMaxNoveltyA"),
                    F.col("maxNoveltyB").alias("shuffledMaxNoveltyB"),
                    "diseaseId",
                    "targetId",
                    "datasourceIdA",
                    "datasourceIdB",
                ),
                ["diseaseId", "targetId", "datasourceIdA", "datasourceIdB"],
                "left",
            )
            .withColumn("realYearA-realYearB", F.col("realYearA") - F.col("realYearB"))
            .withColumn(
                "shuffledYearA-shuffledYearB",
                F.col("shuffledYearA") - F.col("shuffledYearB"),
            )
        ).write.parquet(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/pairs".format(
                evidenceLink
            )
        )

# statistics
if 1:
    for evidenceLink in [
        "direct",
        # "indirect",
    ]:
        data = spark.read.parquet(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/pairs".format(
                evidenceLink
            )
        ).filter(
            ~F.col("datasourceIdA").isin(
                ["crispr_screen", "crispr", "gene_burden", "slapenrich", "sysbio"]
            )
            & ~F.col("datasourceIdB").isin(
                ["crispr_screen", "crispr", "slapenrich", "sysbio", "gene_burden"]
            )
        )

        results_regression_r2_real = []
        results_regression_pval_real = []
        results_ttest_t_real = []
        results_ttest_pval_real = []

        results_regression_r2_random = []
        results_regression_pval_random = []

        for datasourceIdA in (
            data.select("datasourceIdA").distinct().toPandas().datasourceIdA
        ):
            for datasourceIdB in (
                data.select("datasourceIdB").distinct().toPandas().datasourceIdB
            ):
                print(datasourceIdA, datasourceIdB)
                values = data.filter(
                    (F.col("datasourceIdA") == datasourceIdA)
                    & (F.col("datasourceIdB") == datasourceIdB)
                )

                # linear regression for real pairs
                x = values.select("realYearA").toPandas().realYearA.values
                y = values.select("realYearB").toPandas().realYearB.values
                if (len(x) > 1) and (len(y) > 1):
                    rvalue, pvalue = stats.pearsonr(x, y)
                else:
                    rvalue, pvalue = np.nan, np.nan
                results_regression_r2_real.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "r_real": rvalue,
                    }
                )
                results_regression_pval_real.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "p_real": pvalue,
                    }
                )

                # t-test for real pairs
                x = (
                    values.filter(F.col("realYearA-realYearB").isNotNull())
                    .select("realYearA-realYearB")
                    .toPandas()["realYearA-realYearB"]
                    .values
                )
                y = (
                    values.filter(F.col("shuffledYearA-shuffledYearB").isNotNull())
                    .select("shuffledYearA-shuffledYearB")
                    .toPandas()["shuffledYearA-shuffledYearB"]
                    .values
                )
                test = stats.ttest_ind(x, y)
                statistic, pvalue = test.statistic, test.pvalue
                results_ttest_t_real.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "t_real": statistic,
                    }
                )
                results_ttest_pval_real.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "p_real": pvalue,
                    }
                )

                # linear regression for random pairs
                x = (
                    values.filter(~F.col("shuffledYearA").isNull())
                    .select("shuffledYearA")
                    .toPandas()
                    .shuffledYearA.values
                )
                y = (
                    values.filter(~F.col("shuffledYearB").isNull())
                    .select("shuffledYearB")
                    .toPandas()
                    .shuffledYearB.values
                )

                if (len(x) > 1) and (len(y) > 1):
                    rvalue, pvalue = stats.pearsonr(x, y)
                else:
                    rvalue, pvalue = np.nan, np.nan
                results_regression_r2_random.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "r_random": rvalue,
                    }
                )
                results_regression_pval_random.append(
                    {
                        "datasourceIdA": datasourceIdA,
                        "datasourceIdB": datasourceIdB,
                        "p_random": pvalue,
                    }
                )

        results_regression_r2_real = pd.DataFrame(results_regression_r2_real)
        results_regression_r2_real = results_regression_r2_real.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="r_real"
        ).fillna(0)

        results_regression_pval_real = pd.DataFrame(results_regression_pval_real)
        results_regression_pval_real = results_regression_pval_real.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="p_real"
        ).fillna(0)

        results_ttest_t_real = pd.DataFrame(results_ttest_t_real)
        results_ttest_t_real = results_ttest_t_real.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="t_real"
        ).fillna(0)

        results_ttest_pval_real = pd.DataFrame(results_ttest_pval_real)
        results_ttest_pval_real = results_ttest_pval_real.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="p_real"
        ).fillna(0)

        results_regression_r2_random = pd.DataFrame(results_regression_r2_random)
        results_regression_r2_random = results_regression_r2_random.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="r_random"
        ).fillna(0)

        results_regression_pval_random = pd.DataFrame(results_regression_pval_random)
        results_regression_pval_random = results_regression_pval_random.pivot(
            index="datasourceIdA", columns="datasourceIdB", values="p_random"
        ).fillna(0)

        results_regression_r2_real.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/regression_r_real.csv".format(
                evidenceLink
            ),
            sep="\t",
        )
        results_regression_pval_real.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/regression_p_real.tsv".format(
                evidenceLink
            ),
            sep="\t",
        )
        results_ttest_t_real.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/ttest_t_real.tsv".format(
                evidenceLink
            ),
            sep="\t",
        )
        results_ttest_pval_real.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/ttest_p_real.tsv".format(
                evidenceLink
            ),
            sep="\t",
        )

        results_regression_r2_random.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/regression_r_random.tsv".format(
                evidenceLink
            ),
            sep="\t",
        )
        results_regression_pval_random.to_csv(
            "gs://ot-team/cfalaguera/noveltyCorrelation/evidenceLink={}/regression_p_random.tsv".format(
                evidenceLink
            ),
            sep="\t",
        )
