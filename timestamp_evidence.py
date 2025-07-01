#!/usr/bin/env python

__author__ = "Cote Falaguera (mjfalagueramata@gmail.com)"
__date__ = "02 Jun 2025"

import os
import time
from datetime import timedelta
from pyspark.sql import functions as F
from pyspark.sql import types as T
from pyspark.sql import SparkSession
import inspect
import subprocess


"""
timestamp_evidence.py: Annotate Open Targets Platform evidence with their timestamp.
Run this script in Google Cloud after running parse_timestamp.py locally.

Useful GitHub links:
- https://github.com/opentargets/timeseries
- https://github.com/opentargets/issues/issues/2739
"""


# Setup Google Cloud machine and sumbit job:
# gcloud dataproc clusters create cf-timeseries --region europe-west1 --zone europe-west1-d --single-node --master-machine-type n2-highmem-128 --master-boot-disk-size 500 --image-version 2.0-debian10 --project open-targets-eu-dev
# gcloud dataproc jobs submit pyspark timestamp_evidence.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"

# Paths
timestamp_version = "19.05.2025"
timestamp_path = "gs://ot-team/cfalaguera/timestamp/{}/".format(timestamp_version)
ot_platform_version = "25.03"
evidence_path = "gs://open-targets-data-releases/{}/output/evidence".format(
    ot_platform_version
)
credible_set_path = "gs://open-targets-data-releases/{}/output/credible_set".format(
    ot_platform_version
)
study_path = "gs://open-targets-data-releases/{}/output/study".format(
    ot_platform_version
)
evidence_dated_tmp1_path = (
    "gs://ot-team/cfalaguera/evidence_dated_tmp1/"  # source-specific timestamps
)
evidence_dated_tmp2_path = "gs://ot-team/cfalaguera/evidence_dated_tmp2/"  # source-specific timestamps + PubMed timestamps
evidence_dated_path = "gs://ot-team/cfalaguera/{}/evidence_dated/".format(
    ot_platform_version
)
clinvar_timestamp = "gs://ot-team/dsuveges/ClinVar_dated_release.2025.05"  # provisional solution kindly provided by Daniel Suveges

# Establish spark connection
spark = SparkSession.builder.getOrCreate()


# Functions to add timestamps to the evidence files
def date_slapenrich():

    (
        spark.read.parquet(evidence_path + "/sourceId=slapenrich")
        # fill publicationYear and ensure schema conservation
        .fillna({"publicationYear": 2018})
        # fill literature and ensure schema conservation
        .withColumn(
            "literature",
            F.coalesce(
                F.col("literature"), F.array(F.lit("22223"))
            ).cast(  # https://saezlab.github.io/SLAPenrich/, paper publishing the model
                "array<string>"
            ),
        ).write.parquet(evidence_dated_tmp1_path + "/sourceId=slapenrich")
    )


def date_expression_atlas():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    # 3 experiments have forbidden access: E-MTAB-9045, E-MTAB-8887, E-MTAB-9950. They map to the undated pieces of evidence

    source = "expression_atlas"
    (
        spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
        .join(
            spark.read.option("header", "true")
            .csv("{}/{}/evidence2year.csv".format(timestamp_path, source))
            .select(
                "studyId",
                F.col("year").cast(T.IntegerType()).alias("studyYear"),
                F.lit("expression_atlas").alias("datasourceId"),
            ),
            ["studyId", "datasourceId"],
            "left",
        )
        .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
    )


def date_cancer_genome_interpreter():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))

    if 0:
        # do not run this because getYearFromCancerGenomeInterpreter() function messes around with sources field and results in more evidence than originally,
        # until this is fixed, just map pmids in literature to pubmedYear
        (
            spark.read.parquet(evidence_path + "/sourceId=cancer_biomarkers")
            .join(
                spark.read.option("header", "true")
                .csv(
                    timestamp_path
                    + "/sourceId=cancer_genome_interpreter/evidence2year.csv"
                )
                .withColumn("curationYear", F.col("year").cast(T.IntegerType())),
                on=[
                    "targetFromSourceId",
                    "diseaseFromSource",
                    "datasourceId",
                    "drugFromSource",
                    "biomarkerName",
                    "confidence",
                ],
                how="left",
            )
            .write.parquet(evidence_dated_tmp1_path + "/sourceId=cancer_biomarkers")
        )


def date_clingen():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    source = "clingen"
    (
        spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
        .withColumn(
            "url", F.explode_outer("urls.url")
        )  # explode results in only one url per urls
        .join(
            spark.read.option("header", "true")
            .csv("{}/{}/evidence2year.csv".format(timestamp_path, source))
            .withColumn("curationYear", F.col("year").cast(T.IntegerType())),
            on=[
                "targetFromSourceId",
                "diseaseFromSourceId",
                "url",  # never null
                "confidence",
                "studyId",
                "datasourceId",
            ],
            how="left",
        )
        .drop("url")
        .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
    )


def date_gene2phenotype():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    source = "gene2phenotype"
    (
        spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
        .withColumn(
            "tmp", F.explode_outer("allelicRequirements")
        )  # explode results in only one allelicRequirement per allelicRequirements
        .fillna("0", subset=["tmp"])
        .join(
            spark.read.option("header", "true")
            .csv("{}/{}/evidence2year.csv".format(timestamp_path, source))
            .withColumnRenamed("allelicRequirements", "tmp")
            .fillna("0", subset=["tmp"])
            .withColumn("curationYear", F.col("year").cast(T.IntegerType())),
            on=[
                "targetFromSourceId",
                "diseaseFromSource",
                "tmp",
                "confidence",
                "studyId",  # is never null
                "datasourceId",
            ],
            how="left",
        )
        .drop("tmp")
        .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
    )


def date_clinvar():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    for source in ["eva", "eva_somatic"]:
        (
            spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
            .join(
                spark.read.parquet(clinvar_timestamp).select(
                    F.col("rcvAccession").alias("studyId"),
                    F.split("firstSubmissionDate", "-")
                    .getItem(0)
                    .cast(T.IntegerType())
                    .alias("curationYear"),
                ),
                "studyId",  # is never null
                "left",
            )
            .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
        )


def date_gwas_credible_sets():
    # in 25.03 release, publicationYear and literature fields are empty in gwas_credible_sets evidence (they used to be filled in former ot_genetics_portal evidence).
    # to track these new evidence back to the reference literature source, Jack Ge suggested the following solution:
    print("running {} ...".format(inspect.currentframe().f_code.co_name))

    source = "gwas_credible_sets"
    (
        spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
        .drop(
            "studyId"
        )  # in 25.03 studyId is empty for this sourceId=gwas_credible_sets
        .join(
            spark.read.parquet(credible_set_path).select("studyId", "studyLocusId"),
            on="studyLocusId",
            how="left",
        )
        .join(
            spark.read.parquet(study_path).select(
                "studyId",
                F.split(F.col("publicationDate"), "-")
                .getItem(0)
                .cast(T.IntegerType())
                .alias("studyYear"),
            ),
            on="studyId",
            how="left",
        )
        .drop("studyId")
        .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
    )


def date_chembl():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))

    source = "chembl"
    (
        spark.read.parquet("{}/sourceId={}".format(evidence_path, source))
        .withColumn(
            "studyYear",
            F.split(F.col("studyStartDate"), "-").getItem(0).cast(T.IntegerType()),
        )
        .write.parquet("{}/sourceId={}".format(evidence_dated_tmp1_path, source))
    )


def date_pubmed():

    # always process this source the last one
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    source = "pubmed"

    data = spark.read.option("mergeSchema", "true").parquet(
        evidence_dated_tmp1_path
    )  # When reading a directory of Parquet files with potentially different schemas always use mergeSchema option, if not I loose curationYear column cause not all sourceIds have it

    data.unionByName(
        spark.read.parquet(
            evidence_path
        ).join(  # join original evidence not processed with timestamps
            data.select("id"), on="id", how="anti"
        ),
        allowMissingColumns=True,
    ).withColumn("pmid", F.explode_outer("literature")).join(
        spark.read.option("header", "true")
        .csv("{}/{}/evidence2year.csv".format(timestamp_path, source))
        .select(
            F.col("year").cast(T.IntegerType()).alias("pubmedYear"),
            F.col("pmid")
            .cast("string")
            .alias(
                "pmid"
            ),  # pmids within 'literature' field on OTP evidence files are strings, and we will join them
        ),
        on="pmid",
        how="left",
    ).write.partitionBy(
        "sourceId"
    ).parquet(
        evidence_dated_tmp2_path
    )


def select_best_date():
    """
    Choose best timestamp found for each evidence: studyYear > curationYear > publicationYear > pubmedYear
    """

    print("running {} ...".format(inspect.currentframe().f_code.co_name))

    data = spark.read.parquet(evidence_dated_tmp2_path)

    (
        data.unionByName(
            spark.read.parquet(evidence_path)
            .join(  # join original evidence not processed with timestamps
                data.select("id"), on="id", how="anti"
            )
            .withColumn("pmid", F.explode_outer("literature")),
            allowMissingColumns=True,
        )
        .withColumn(
            "year",
            F.coalesce("studyYear", "curationYear", "publicationYear", "pubmedYear"),
        )  # coalesce() checks columns in the specified order and returns the first non-null value
        .groupBy(
            "id",
            "targetId",
            "diseaseId",
            "datasourceId",
            "score",
            "drugId",
            "clinicalPhase",
            "studyLocusId",
            "urls",
        )
        .agg(F.collect_set("pmid").alias("literature"), F.min("year").alias("year"))
        .withColumn("sourceId", F.col("datasourceId"))
        .write.partitionBy("sourceId")
        .parquet(evidence_dated_path)
    )

    # remove intial tmp files
    command = ["gsutil", "rm", "-r", evidence_dated_tmp1_path, evidence_dated_tmp2_path]
    result = subprocess.run(command, capture_output=True, text=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    print("Return code:", result.returncode)


def count_dated_evidence():
    """
    Count evidence by source with or without year found.

    Returns:
        Dataframe with evidence mapped to their publication year. Columns:
        - datasourceId
        - evidence
        - dated_evidence
    """

    (
        spark.read.parquet(evidence_dated_path)
        .groupBy("datasourceId")
        .agg(F.size(F.collect_set("id")).alias("evidence"))
        .join(
            spark.read.parquet(evidence_dated_path)
            .filter(F.col("year").isNotNull())
            .groupBy("datasourceId")
            .agg(F.size(F.collect_set("id")).alias("dated_evidence")),
            ["datasourceId"],
            "left",
        )
        .show(50, truncate=False)
    )


# Run functions
start_time = time.perf_counter()

date_slapenrich()
date_expression_atlas()
date_cancer_genome_interpreter()
date_clingen()
date_gene2phenotype()
date_clinvar()
date_gwas_credible_sets()
date_chembl()
date_pubmed()  # always run this the last one

select_best_date()

count_dated_evidence()

end_time = time.perf_counter()
elapsed_seconds = end_time - start_time
elapsed_td = timedelta(seconds=elapsed_seconds)
days = elapsed_td.days
hours, remainder = divmod(elapsed_td.seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"\nElapsed time: {days:02d}-{hours:02d}:{minutes:02d}:{seconds:02d}")

# Timestamps for the rest of sources are extracted from PubMed PMIDs when available:
# europepmc - publicationYear is never null! :)
# sysbio - 100% covered by pmid -> map to PubMed years
# crispr_brain - some are covered by pmid -> map to PubMed years
# eva and eva_somatic provisional solution: gs://ot-team/dsuveges/ClinVar_dated_release.2023.07 (eva_somatic 100% covered, eva: 2739522 out of 2742321)
# intogene - releases since 2013 no mention to publication date
# orphanet - unsuccessfully tried to find dates in associations file https://www.orphadata.com/genes/
# progeny - no date - @irene: What we call Progeny is the application of this algorithm to gene expression data from TCGA. The data is here https://storage.googleapis.com/otar000-evidence_input/PROGENy/data_files/progeny_normalVStumor_opentargets.txt
# genomics_england - @irene's link doesn't have version year: https://storage.googleapis.com/otar000-evidence_input/GenomicsEngland/data_files/All_genes_20220804-1350_green_amber_red_public_v1plus_no_superpanels.tsv
# cancer_gene_census - PubMed PMID
# reactome - PubMed PMID
# impc - no mention to date in files
# gene_bass already provides literature pmid: 36778668 from https://app.genebass.org/terms
# pubmed and clinvar too large to run them locally
# project score evidence where mapped in 23.06 to pmid 30971826 (2019, Behan) and in 25.03 are mapped to pmid 38215750 (2024, Pacini). I can't track back which evidence come from each study...
