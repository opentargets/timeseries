#!/usr/bin/env python

__author__ = "Cote Falaguera (mjfalagueramata@gmail.com)"
__date__ = "02 Jul 2025"

"""
timeseries.py: Assess the evolution over time of evidence supporting target-disease associations in the Open Targets Platform.

Useful GitHub links:
- https://github.com/opentargets/timeseries
- https://github.com/opentargets/issues/issues/2739
"""
import datetime
import time
from datetime import timedelta
from typing import Any
import yaml

from pyspark.sql import functions as f
from pyspark.sql import SparkSession, Window, DataFrame, Column
from pyspark.sql import types as t



def get_weight_for_datasource(
        data_sources: dict[str, Any]
) -> list[tuple[str, float]]:
    """Get list of data sources' weights for overall score.

    Args:
        data_sources (dict[str, Any]): 
    """
    return [(datasource["id"], datasource["weight"]) for datasource in data_sources]

# Novelty assessment functions
def get_indirect_evidence(
        evidence_df: DataFrame,
        disease_df: DataFrame,
):
    """
    Propagate dated evidence accross disease ontology.
    """
    (
        evidence_df
        .join(
            disease_df.select(
                f.col("id").alias("diseaseId"),
                f.explode(
                    f.array_union(f.array(f.col("id")), f.col("ancestors"))
                ).alias("specificDiseaseId"),
            ),
            "diseaseId",
            "inner",
        )  # add ancestors
        .drop("diseaseId")  # drop descendants
        .withColumnRenamed("specificDiseaseId", "diseaseId")
    )

def evidence_date_to_year(evidence_date: Column, last_year: int) -> Column:
    """Convert evidence date to year when it is not available, use the last year +1
    
    Args:
        eivdence_date (Column): spark column with the evidence date
        last_year (int): last year considered in analysis
    
    Returns:
        Column: year, which should never be null.
    """
    return (
        f.when(
            evidence_date.isNotNull(), 
            f.split(evidence_date, '-')[0].cast(t.IntegerType())
        ).otherwise(
            f.lit(last_year)
        )
    )

def get_harmonic_sum(collected_scores: Column, max_value: int) -> Column:
    """Calculating harmonic sum on an array of score values
    
    Process:
        1. Sorting collected scores
        2. Get top 50 scores
        3. get weighted indices 
        4. compute harmonic sum of scores
        5. Normalise harmonic sum

    Args:
        collected_scores (Column): list of floats in a spark columns
    
    Returns:
        Column: float normalised harmonic sum of the input array
    """
    sorted_filtered_scores = f.slice(
        f.reverse(f.array_sort(collected_scores)),
        1, 50
    ).alias('scores')

    indices = f.sequence(
        f.lit(1), f.size(sorted_filtered_scores)
    ).alias('indices')
    return indices
    # weighted_scores = f.transform(
    #     f.arrays_zip(sorted_filtered_scores, indices),
    #     lambda pair: pair.scores / f.pow(pair.indices, f.lit(2))
    # )

    # return f.sum(weighted_scores) / max_value

def get_association_score_by_datasource_dated(
    evidence_df: DataFrame,
    last_year: str,
    first_year: str,
    max_score: float,
    excluded_datasources:list[str]= [],
    excluded_diseases:list[str]=[ ],
    excluded_targets:list[str]= [],
) -> DataFrame:
    """
    Recalculate association scores by datasource over the years. Estimated running time: 2h.

    Args:
        evidence_df (DataFrame): Input Spark DataFrame containing evidence records.
        last_year (str): The last year to include in the analysis.
        first_year (str): The first year to include in the analysis.
        max_score (float): Maximum score value for normalization.
        excluded_datasources (list[str], optional): List of datasource IDs to exclude.
        excluded_diseases (list[str], optional): List of disease IDs to exclude.
        excluded_targets (list[str], optional): List of target IDs to exclude.

    Returns:
        Dataframe with association scores by datasource over the years.
    """

    # Exclude data_source if provided:
    if len(excluded_datasources) > 0:
        evidence_df = evidence_df.filter(~f.col("datasourceId").isin(excluded_datasources))
    
    # Exclude diseases if provided:
    if len(excluded_diseases) > 0:
        evidence_df = evidence_df.filter(~f.col("diseaseId").isin(excluded_diseases))

    # Exclude targets if provided:
    if len(excluded_targets) > 0:
        evidence_df = evidence_df.filter(~f.col("targetId").isin(excluded_targets))


    # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
    evidence_df = (
        evidence_df.withColumn(
            "year", 
            evidence_date_to_year(f.col('evidenceDate'), last_year)
        )
    )

    # get all the combinations of datasourceId vs years in the range between the first_year and the last_year set
    sourceVSyear = (
        # unique sources
        evidence_df.select("datasourceId")
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
    evidence_df = (
        sourceVSyear.join(
            evidence_df.select("diseaseId", "targetId", "datasourceId").distinct(),
            "datasourceId",
            "left",
        )
        # disease - target - datasource - year - score
        .join(
            evidence_df,
            ["diseaseId", "targetId", "datasourceId", "year"],
            "left",
        )
    )

    # prepare partition: all evidence accumulated for each disease-target-datasource-year triplet until the given year:
    partition1 = (
        Window.partitionBy("diseaseId", "targetId", "datasourceId")
        .orderBy("year")
        .rangeBetween(Window.unboundedPreceding, 0)
    )

    # recalculate harmonic sum scores over the years considering the evidence accumulated:
    return (
        evidence_df.groupBy("diseaseId", "targetId", "datasourceId", "year")
        # collect list of scores FOR each year:
        .agg(f.collect_list("score").alias("cum_scores"))
        # collect scores UNTIL each year:
        .withColumn(
            "cum_scores", f.flatten(f.collect_list("cum_scores").over(partition1))
        )
        .select(
            'targetId',
            'diseaseId',
            'datasourceId',
            'cum_scores',
            'year'
            get_harmonic_sum(f.col('cum_scores')).alias('harmonicSum')
        )
        .show()
    )
        # # remove NaNs from the cumulative scores array
        # .withColumn(
        #     "scores_no_nan", f.expr("filter(cum_scores, x -> NOT isnan(x))")
        # )
        # sort descending and take top 50 scores
        .withColumn("scores_sorted", f.reverse(f.array_sort("scores_no_nan")))
        .withColumn("top50_scores", f.expr("slice(scores_sorted, 1, 50)"))
        # generate indices (1-based) for the top 50 scores
        .withColumn("idx", f.sequence(f.lit(1), f.size("top50_scores")))
        # zip scores and indices, and divide each score by idx^2
        .withColumn(
            "weighted_scores",
            f.expr(
                "transform(arrays_zip(top50_scores, idx), x -> x.top50_scores / pow(x.idx, 2))"
            ),
        )
        # sum the weighted scores
        .withColumn(
            "harmonic_sum",
            f.expr("aggregate(weighted_scores, 0D, (acc, x) -> acc + x)"),
        )
        # normalize by max_score
        .withColumn(
            "harmonic_score", 
            f.col("harmonic_sum") / f.lit(max_score)
        )
        # select
        .select(
            "targetId",
            "diseaseId",
            "year",
            "datasourceId",
            f.col("harmonic_score").alias("score"),
        )
        # recover non-dated evidence
        .withColumn(
            "year",
            f.when(f.col("year") == last_year + 1, None).otherwise(f.col("year")),
        )
    )


# def get_association_novelty_by_datasource_dated(
#     evidenceLink,
#     scale:float,
#     shift:float,
#     window:float,
#     diseaseId=[],
#     targetId=[],
#     excludeDatasource=[],
# ):
#     """
#     Calculate novelty of association scores by datasource over the years.

#     Args:
#         evidenceLink (str): 'direct' or 'indirect'
#         scale (float):      logistic growth rate or steepness of the novelty curve
#         shift (float):      x-axis value of the sigmoid's point of the novelty curve
#         window (float):     range of years after the peak year to apply decay function to

#     Returns:
#         Dataframe with novelty of association scores by datasource over the years. Columns:
#         - diseaseId
#         - targetId
#         - datasourceId
#         - year
#         - score
#         - novelty
#     """

#     # get dated score
#     if evidenceLink == "direct":
#         f = association_by_datasource_dated_file
#     elif evidenceLink == "indirect":
#         f = association_by_datasource_dated_indirect_file
#     data = spark.read.parquet(f)

#     # exclude data_source
#     if len(excludeDatasource):
#         data = data.filter(~f.col("datasourceId").isin(excludeDatasource))

#     if len(diseaseId):
#         data = data.filter(~f.col("diseaseId").isin(diseaseId))

#     if len(targetId):
#         data = data.filter(~f.col("targetId").isin(targetId))

#     # prepare partition: disease-target-datasource triplet ordered by increasing year
#     partition1 = Window.partitionBy("diseaseId", "targetId", "datasourceId").orderBy(
#         "year"
#     )

#     # calculate novelty for disease-target-datasource over the years
#     data = (
#         data
#         # fill non-dated evidence with last_year+1 to avoid loosing them (we'll reset them later)
#         .fillna(last_year + 1, subset=["year"])
#         # fill NaN score with 0 for later novelty calculation
#         .fillna(0, subset=["score"])
#         # for each target-disease-datasource get peaks of score shift (when current year score minus previous year score > 0)
#         .select(
#             "diseaseId",
#             "targetId",
#             "datasourceId",
#             f.col("year").alias("peakYear"),
#             (f.col("score") - f.lag("score", offset=1).over(partition1)).alias("peak"),
#         )
#         # filter peaks
#         .filter(f.col("peak") > 0)
#         # for each peak year, get the range of years within a window, e.g. peakYear=1991 -> range of peakYear + window=3: 1991, 1992, 1993
#         .select(
#             "*",
#             f.posexplode(  # Returns a new row for each element with position in the given array or map. Uses the default column name pos for position, and col for elements in the array and key and value for elements in the map unless specified otherwise.
#                 f.sequence(  # Generate a sequence of integers from start to stop, incrementing by step. If step is not set, incrementing by 1 if start is less than or equal to stop, otherwise -1.
#                     f.col("peakYear"),
#                     f.col("peakYear") + window,
#                 )
#             ).alias(
#                 "year-peakYear", "year"
#             ),
#         )
#         # for each peak, calculate the novelty value at the different years within the window (novelty = peakScore/(1+exp^(scale*(year-peakYear-shift))) -> logistic function)
#         # and select max. novelty value found for each year
#         .groupBy("diseaseId", "targetId", "datasourceId", "year")
#         .agg(
#             f.round(
#                 f.max(
#                     f.col("peak")
#                     / (1 + f.exp(scale * (f.col("year-peakYear") - shift))),
#                 ),
#                 3,
#             ).alias("novelty")
#         )
#         # add max. novelty values to original disease-target-datasource-year dataframe
#         .join(
#             data,
#             ["diseaseId", "targetId", "datasourceId", "year"],
#             "right",
#         )
#         # reset non-dated evidence to null year
#         .replace(last_year + 1, None, subset=["year"])
#         # set novelty=0 when novelty=null
#         .fillna(0, subset=["novelty"])
#         .withColumn("sourceId", f.col("datasourceId"))
#         .write.mode("overwrite")
#         .partitionBy("sourceId")
#         .parquet(f)
#     )


# def get_association_score_by_overall_dated(
#     evidenceLink="indirect",
#     excludeDatasource=[],
# ):
#     """
#     Recalculate overall association scores over the years.

#     Args:
#         evidenceLink (str):         'direct' or 'indirect'
#         excludeDatasource (list):   data_source to exclude, e.g. "chembl"

#     Returns:
#         Dataframe with overall association scores over the years. Columns:
#         - diseaseId
#         - targetId
#         - year
#         - score
#     """

#     # get dated score
#     if evidenceLink == "direct":
#         data = spark.read.parquet(association_by_datasource_dated_file)
#         f = association_by_overall_dated_file
#     elif evidenceLink == "indirect":
#         data = spark.read.parquet(association_by_datasource_dated_indirect_file)
#         f = association_by_overall_dated_indirect_file

#     if os.path.exists(f):
#         pass

#     else:

#         # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
#         data = data.fillna(last_year + 1, subset=["year"])

#         # get data_source' weights in overall score formula
#         weights = spark.createDataFrame(
#             data=[
#                 [datasourceId, str(weight)]
#                 for datasourceId, weight in get_weight_for_datasource()
#             ],
#             schema=["datasourceId", "weight"],
#         )

#         # exclude data_source
#         if len(excludeDatasource):
#             data = data.filter(~f.col("datasourceId").isin(excludeDatasource))

#         # recalculate harmonic sum scores over the years considering the evidence accumulated
#         data = (
#             # add data_source' weights
#             data.join(weights, "datasourceId", "left")
#             # weight source-specific scores
#             .withColumn("score", f.col("score") * f.col("weight"))
#             # list of weighted scores
#             .groupBy("diseaseId", "targetId", "year")
#             # collect list of scores FOR each year
#             .agg(f.collect_list("score").alias("cum_scores"))
#             # remove NaNs from the cumulative scores array
#             .withColumn(
#                 "scores_no_nan", f.expr("filter(cum_scores, x -> NOT isnan(x))")
#             )
#             # sort descending and take top 50 scores
#             .withColumn("scores_sorted", f.reverse(f.array_sort("scores_no_nan")))
#             .withColumn("top50_scores", f.expr("slice(scores_sorted, 1, 50)"))
#             # generate indices (1-based) for the top 50 scores
#             .withColumn("idx", f.sequence(f.lit(1), f.size("top50_scores")))
#             # zip scores and indices, and divide each score by idx^2
#             .withColumn(
#                 "weighted_scores",
#                 f.expr(
#                     "transform(arrays_zip(top50_scores, idx), x -> x.top50_scores / pow(x.idx, 2))"
#                 ),
#             )
#             # sum the weighted scores
#             .withColumn(
#                 "harmonic_sum",
#                 f.expr("aggregate(weighted_scores, 0D, (acc, x) -> acc + x)"),
#             )
#             # normalize by max_score
#             .withColumn("harmonic_score", f.col("harmonic_sum") / f.lit(max_score))
#             # select
#             .select(
#                 "targetId",
#                 "diseaseId",
#                 "year",
#                 f.col("harmonic_score").alias("score"),
#             )
#             # recover non-dated evidence
#             .withColumn(
#                 "year",
#                 f.when(f.col("year") == last_year + 1, None).otherwise(f.col("year")),
#             )
#             .write.parquet(f)
#         )


# def get_association_novelty_by_overall_dated(
#     evidenceLink="indirect",
#     scale=novelty_scale,
#     shift=novelty_shift,
#     window=novelty_window,
#     excludeDatasource=[],
# ):
#     """
#     Calculate overall novelty of association scores over the years.

#     Args:
#         evidenceLink (str):         'direct'/'indirect'
#         scale (float):              logistic growth rate or steepness of the novelty curve
#         shift (float):              x value of the sigmoid's point of the novelty curve
#         window (float):             range of years after the peak year to apply decay function to
#         excludeDatasource (list):   data_source to exclude, e.g. "chembl"

#     Returns:
#         Dataframe with overall novelty of association scores over the years. Columns:
#         - diseaseId
#         - targetId
#         - year
#         - score
#         - novelty
#     """

#     # get dated score
#     if evidenceLink == "direct":
#         f = association_by_overall_dated_file
#     elif evidenceLink == "indirect":
#         f = association_by_overall_dated_indirect_file
#     data = spark.read.parquet(
#         f
#     ).persist()  # persist since we will call this dataframe more than once

#     # fill non-dated evidence with last_year + 1 to avoid loosing them (we'll reset them later)
#     data = data.fillna(last_year + 1, subset=["year"])

#     print(data.show())

#     # calculate novelty for disease-target over the years
#     (
#         data
#         # fill NaN score with 0 for later novelty calculation
#         .fillna(0, subset=["score"])
#         # for each target-disease get peaks of score shift (when current year score minus previous year score > 0)
#         .select(
#             "diseaseId",
#             "targetId",
#             f.col("year").alias("peakYear"),
#             #  partition: disease-target pair ordered by increasing year
#             (
#                 f.col("score")
#                 - f.lag("score", offset=1).over(
#                     Window.partitionBy("diseaseId", "targetId").orderBy("year")
#                 )
#             ).alias("peak"),
#         )
#         # filter peaks
#         .filter(f.col("peak") > 0)
#         # for each peak year, get the range of years within a window, e.g. peakYear=1991 -> range of peakYear + window=3: 1991, 1992, 1993
#         .select(
#             "*",
#             f.posexplode(
#                 f.sequence(
#                     f.col("peakYear"),
#                     f.col("peakYear") + window,
#                 )
#             ).alias("year-peakYear", "year"),
#         )
#         # for each peak, calculate the novelty value at the different years within the window (novelty = peakScore/(1+exp^(scale*(year-peakYear-shift))) -> logistic function)
#         # and select max. novelty value found for each year
#         .groupBy("diseaseId", "targetId", "year")
#         .agg(
#             f.round(
#                 f.max(
#                     f.col("peak")
#                     / (1 + f.exp(scale * (f.col("year-peakYear") - shift))),
#                 ),
#                 3,
#             ).alias("novelty")
#         )
#         # add max. novelty values to original dataframe disease-target-year
#         .join(
#             data,
#             ["diseaseId", "targetId", "year"],
#             "right",
#         )
#         # reset non-dated evidence to null year
#         .replace(last_year + 1, None, subset=["year"])
#         # set novelty=0 when novelty=null
#         .fillna(0, subset=["novelty"])
#         .write.mode("overwrite")
#         .parquet(f)
#     )

#     print(data.show())
#     data.unpersist()

def report_elapsed_time(start_time: datetime) -> None:
    """Print out elapsed time to log the time it took to run the script.
    
    Args:
        start_time (datetime): timestamp when the process started.
    """
    end_time = time.perf_counter()
    elapsed_seconds = end_time - start_time
    elapsed_td = timedelta(seconds=elapsed_seconds)
    days = elapsed_td.days
    hours, remainder = divmod(elapsed_td.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\nElapsed time: {days:02d}-{hours:02d}:{minutes:02d}:{seconds:02d}")

def read_yaml_config(path):
    with open(path, 'r') as f:
        config = yaml.safe_load(f)
    return config

if __name__ == '__main__':

    # Start measuring time:
    start_time = time.perf_counter()

    # Establish spark connection
    spark = (
        SparkSession
        .builder
        .config("spark.driver.memory", "10G")
        .getOrCreate()
    )

    INPUT_PATH = '/Users/dsuveges/project_data/25.09/output/'
    OUTPUT_PATH = '/Users/dsuveges/project_data/25.09/view/'

    # Reading config:
    config: dict = read_yaml_config('timeseries_configuration.yaml')

    # Extract configuration:
    first_year = config.get('first_year')
    novelty_scale = config.get('novelty_scale')
    novelty_shift = config.get('novelty_shift')
    novelty_window = config.get('novelty_window')
    max_score = config.get('max_score')

    # Last year is dynamic:
    last_year = datetime.datetime.today().year

    # Reading necessary files:
    direct_evidence_df = spark.read.parquet(f'{INPUT_PATH}/evidence/sourceId=gwas_credible_sets')
    disease_df = spark.read.parquet(f'{INPUT_PATH}/disease')

    # Explode evidence to all ancestors of the associated disease:
    # indirect_evidence_df = get_indirect_evidence(
    #     direct_evidence_df,
    #     disease_df
    # )

    (
        get_association_score_by_datasource_dated(
            direct_evidence_df, last_year, first_year, max_score
        )
        .write.mode('overwrite')
        .parquet(f'{OUTPUT_PATH}/novelty_by_datasource')
    )

    # get_association_novelty_by_datasource_dated(
    #     evidenceLink="direct",
    # )
    # get_association_score_by_datasource_dated(
    #     evidenceLink="indirect",
    # )
    # get_association_novelty_by_datasource_dated(
    #     evidenceLink="indirect",
    # )
    # get_association_score_by_overall_dated(evidenceLink="direct")
    # get_association_novelty_by_overall_dated(evidenceLink="direct")
    # get_association_score_by_overall_dated(evidenceLink="indirect")
    # get_association_novelty_by_overall_dated(evidenceLink="indirect")

    # Report elapsed time:
    report_elapsed_time(start_time)