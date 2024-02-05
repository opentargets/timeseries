#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "05 Feb 2024"

import xml.etree.ElementTree as ET
import gzip
import os
import ftplib
import pandas as pd
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
import numpy as np

"""
timestampsParser.py: Parse timestamps for OT Platform evidence from original resources.
"""

# useful github links:
# https://github.com/opentargets/issues/issues/2739

# Paths
db_path = "/Users/mariaf/"


def getYearFromGeneBass(version="2022/04/13"):
    # data_path = "/Users/mariaf/GeneBass/{}/".format(version)
    # os.system(
    #    "gsutil cp -r gs://ukbb-exome-public/500k/results/results.mt {}/." # https://app.genebass.org/downloads
    # )
    # hl.read_matrix_table('/Users/mariaf/GeneBass/2022.04.13/results.mt')
    # no field with date info

    # Genebass is the only projectId in datasourceId=gene_burden with no literature annotation
    return {
        "datasourceId": "gene_burden",
        "projectId": "Genebass",
        "literature": [36778668],  # https://app.genebass.org/terms (cite us)
        # "year": 2022,
    }


def getYearFromSlapenrich():
    # all slapenrich evidence lack literature annotation
    return {
        "datasourceId": "slapenrich",
        "literature": [
            29713020
        ],  # https://saezlab.github.io/SLAPenrich/, paper publishing the model
        # "year": 2018,
    }


def getYearFromCRISPRbrain():
    return {
        "datasourceId": "crispr_screen",
        "projectId": "crispr_brain",
        "literature": [
            34031600
        ],  # https://crisprbrain.org/about/ (preprint), https://www.nature.com/articles/s41593-021-00862-0 (publication)
        # "year": 2021,
    }


def getYearFromProjectScore():
    return {
        "datasourceId": "crispr",
        "literature": [
            30971826
        ],  # https://score.depmap.sanger.ac.uk/documentation, Fiona behan supp mat
        # "year": 2019,
    }


def getYearFromPubMedBaseline(version="24.07.2023"):
    # yearly update
    url = "ftp.ncbi.nlm.nih.gov"
    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd("pubmed/baseline")

    pubmed_path = "{}/PubMed/{}/".format(db_path, version)

    tmpFile = pubmed_path + "/baseline/{}"
    outputFile = pubmed_path + "/baseline/pmid2year/{}"

    for inputFile in ftp.nlst("*.xml.gz"):
        print(inputFile)
        if not os.path.exists(outputFile.format(inputFile.replace(".xml.gz", ".csv"))):
            ftp.retrbinary(
                "RETR " + inputFile, open(tmpFile.format(inputFile), "wb").write
            )

            try:
                with gzip.open(tmpFile.format(inputFile), "rt") as fht, open(
                    outputFile.format(inputFile.replace(".xml.gz", ".csv")), "w"
                ) as fho:
                    fho.write("pmid,year\n")
                    tree = ET.parse(fht)
                    root = tree.getroot()
                    for i in root.findall("./PubmedArticle/PubmedData"):
                        for j in i.findall(
                            'History/PubMedPubDate[@PubStatus="pubmed"]/Year'
                        ):
                            year = int(j.text)
                        for j in i.findall('ArticleIdList/ArticleId[@IdType="pubmed"]'):
                            pmid = int(j.text)
                        fho.write("{},{}\n".format(pmid, year))
                    print(fho.name, "generated")

                os.remove(tmpFile.format(inputFile))

            except:
                os.remove(outputFile.format(inputFile.replace(".xml.gz", ".csv")))
                os.remove(tmpFile.format(inputFile))
                print(inputFile, "error")

    ftp.quit()


def getYearFromPubMedUpdatefiles(version="24.07.2023"):
    # daily update
    url = "ftp.ncbi.nlm.nih.gov"
    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd("pubmed/updatefiles")

    pubmed_path = "{}/PubMed/{}/".format(db_path, version)

    tmpFile = pubmed_path + "/updatefiles/{}"
    outputFile = pubmed_path + "/updatefiles/pmid2year/{}"

    for inputFile in ftp.nlst("*.xml.gz"):
        print(inputFile)
        if not os.path.exists(outputFile.format(inputFile.replace(".xml.gz", ".csv"))):
            ftp.retrbinary(
                "RETR " + inputFile, open(tmpFile.format(inputFile), "wb").write
            )

            try:
                with gzip.open(tmpFile.format(inputFile), "rt") as fht, open(
                    outputFile.format(inputFile.replace(".xml.gz", ".csv")), "w"
                ) as fho:
                    fho.write("pmid,year\n")
                    tree = ET.parse(fht)
                    root = tree.getroot()
                    for i in root.findall("./PubmedArticle/PubmedData"):
                        for j in i.findall(
                            'History/PubMedPubDate[@PubStatus="pubmed"]/Year'
                        ):
                            year = int(j.text)
                        for j in i.findall('ArticleIdList/ArticleId[@IdType="pubmed"]'):
                            pmid = int(j.text)
                        fho.write("{},{}\n".format(pmid, year))
                    print(fho.name, "generated")

                os.remove(tmpFile.format(inputFile))

            except:
                os.remove(outputFile.format(inputFile.replace(".xml.gz", ".csv")))
                os.remove(tmpFile.format(inputFile))
                print(inputFile, "error")

    ftp.quit()


# not found in https://github.com/opentargets/evidence_datasource_parsers/tree/master/modules
def getYearFromExpressionAtlas(version=2023.06):
    # for the experiments annotated with a pmid in evidence file, "Public Release Date" matches the pmid publication year, and is not empty for those lacking pmid
    url = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd("pub/databases/microarray/data/atlas/experiments")

    expressionatlas_path = "{}/ExpressionAtlas/{}/".format(
        db_path, version
    )  # Pedro Madrigal, contact person, http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments

    tmpFile = expressionatlas_path + "/experiments/{}"
    outputFile = expressionatlas_path + "/experiment2year.csv"

    with open(outputFile, "w") as fho:
        fho.write("studyId,year,datasourceId\n")
        for inputFile in ftp.nlst("E-*"):
            print(inputFile)

            try:
                ftp.retrbinary(
                    "RETR " + "{}/{}.idf.txt".format(inputFile, inputFile),
                    open(tmpFile.format(inputFile), "wb").write,
                )
            # exepcional cases are named <experiment> instead of <experiment>/<experiment>.idf.txt
            except ftplib.error_perm:
                try:
                    ftp.retrbinary(
                        "RETR " + inputFile,
                        open(tmpFile.format(inputFile), "wb").write,
                    )
                # exceptional cases have empty folder
                except ftplib.error_perm:
                    print(inputFile, "error")
                    os.remove(tmpFile.format(inputFile))
                    continue

            try:
                with open(tmpFile.format(inputFile), "r") as fht:
                    studyId = inputFile
                    for line in fht:
                        if "Public Release Date" in line:
                            year = int(
                                line.split("\t")[1].split("-")[0].strip('"')
                            )  # some files use quotes
                            fho.write(
                                "{},{},{}\n".format(studyId, year, "expression_atlas")
                            )
                            break
            # exceptional cases are encoded
            except UnicodeDecodeError:
                with open(tmpFile.format(inputFile), "rb") as fht:
                    studyId = inputFile
                    for line in fht:
                        line = line.decode()
                        if "Public Release Date" in line:
                            year = int(
                                line.split("\t")[1]
                                .split("-")[0]
                                .strip('"')  # exeptional cases use quotes
                            )
                            fho.write("{},{}\n".format(studyId, year))
                            break

            os.remove(tmpFile.format(inputFile))
        print(fho.name, "generated")

    ftp.quit()


def tranformDateForCancerGenomeInterpreter(date):
    if date != date:  # nan
        year = None
    elif "." in date:
        year = int(date.split(".")[-1])
    elif "/" in date:
        year = int("20" + date.split("/")[-1])
    return year


# easy: take Curation date in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/cancerBiomarkers.py
def getYearFromCancerGenomeInterpreter(version="2018.01.17"):
    cancergenomeinterpreter_path = "{}/CancerGenomeInterpreter/{}/".format(
        db_path,
        version,
    )  # https://www.cancergenomeinterpreter.org/2018/data/cgi_biomarkers_latest.zip -> https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/cancerBiomarkers.py

    # datasourceId=cancer_biomarkers
    data = pd.read_csv(
        cancergenomeinterpreter_path + "cgi_biomarkers_20221006.tsv", sep="\t"
    )[
        [
            "Gene",
            "Curation date",
            "Primary Tumor type full name",
            "Drug",
            "Biomarker",
            "Evidence level",
            "Drug full name",
            "Source",
        ]
    ]
    # split and explode drugs
    data["Drug"] = data["Drug"].fillna("")
    data["Drug"] = data.apply(
        lambda row: row["Drug"]
        .strip("[")
        .strip("]")
        .replace(";", ",")
        .split(
            ","
        ),  # sometimes they use semicolon-separated lists ([]) and some others comma-separated strings
        axis=1,
    )
    data = data.explode("Drug")
    data["Drug full name"] = data["Drug full name"].fillna("")
    data["Drug full name"] = data.apply(
        lambda row: row["Drug full name"].strip("[").strip("]").split(";"), axis=1
    )  # in this case, only semicolon-separated strings have been found, which may conain inside commas that we don't want to remove
    data = data.explode("Drug full name")

    # in evidence file, when no Drug is found, Drug full name is used
    data["drugFromSource"] = data.apply(
        lambda row: row["Drug full name"] if row["Drug"] == "" else row["Drug"],
        axis=1,
    )
    # capitalise each letter in the phrase to make it match with evidence table
    data["drugFromSource"] = data.apply(
        lambda row: row["drugFromSource"].title(), axis=1
    )
    # "source" field is a mess, sometimes it can be only a pmid, sometimes it is a semicolon-separated list of url and PMID preceded by "PMID:" tag, or "PMID:  " tag, thus.
    # Thus, I ignore this field
    data = data.rename(
        columns={
            "Gene": "targetFromSourceId",
            "Primary Tumor type full name": "diseaseFromSource",
            "Curation date": "curationDate",
            "Biomarker": "biomarkerName",
            "Evidence level": "confidence",
        }
    ).drop(columns=["Drug full name", "Drug", "Source"])
    data["year"] = data.apply(
        lambda row: tranformDateForCancerGenomeInterpreter(date=row["curationDate"]),
        axis=1,
    )
    data = data.drop(columns=["curationDate"])
    data["datasourceId"] = "cancer_biomarkers"
    data.to_csv(cancergenomeinterpreter_path + "evidence2year.csv", index=False)
    print(cancergenomeinterpreter_path + "evidence2year.csv", "generated")


# easy: use CLASSIFICATION DATE which is already extracted in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/ClinGen.py
def getYearFromClinGen(version="2023.07.25"):
    # https://search.clinicalgenome.org/kb/gene-validity/download

    # 27 very recently curated associations are missing in the raw file but appear in the website

    clingen_path = "{}/ClinGen/{}/".format(db_path, version)
    with open(
        clingen_path
        + "Clingen-Gene-Disease-Summary-{}.csv".format(version.replace(".", "-")),
        "r",
    ) as fhi, open(
        clingen_path
        + "Clingen-Gene-Disease-Summary-{}_wHeader.csv".format(
            version.replace(".", "-")
        ),
        "w",
    ) as fho:
        i = 0
        while i < 4:
            next(fhi)
            i += 1
        for line in fhi:
            print(line)
            fho.write(line)

    data = pd.read_csv(
        clingen_path
        + "Clingen-Gene-Disease-Summary-{}_wHeader.csv".format(
            version.replace(".", "-")
        )
    ).rename(
        columns={
            "GENE SYMBOL": "targetFromSourceId",
            "DISEASE ID (MONDO)": "diseaseFromSourceId",
            "CLASSIFICATION": "confidence",
            "ONLINE REPORT": "url",
            "CLASSIFICATION DATE": "year",
            "GCEP": "studyId",
        }
    )[
        [
            "targetFromSourceId",
            "diseaseFromSourceId",
            "confidence",
            "url",
            "year",
            "studyId",
        ]
    ]
    data["year"] = data.apply(lambda row: row.year.split("-")[0], axis=1)
    data["datasourceId"] = "clingen"
    data.to_csv(clingen_path + "evidence2year.csv", index=False)
    print(clingen_path + "evidence2year.csv")


# easy: use gene disease pair entry date which is already extracted in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/Gene2Phenotype.py
def getYearFromGene2Phenotype(version="28_04_2023"):
    # some of them have pmid annotated, some other don't but all of them have "gene disease pair entry date", shuld we prioritise the second one?
    url = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd("pub/databases/gene2phenotype/{}".format(version))

    gene2phenotype_path = "{}/Gene2Phenotype/{}/".format(
        db_path, version
    )  # http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/28_04_2023/

    tmpFile = gene2phenotype_path + "/{}"
    outputFile = gene2phenotype_path + "/pair2year/{}"

    for inputFile in ftp.nlst("*.csv.gz"):
        print(inputFile)
        if not os.path.exists(outputFile.format(inputFile.replace(".csv.gz", ".csv"))):
            ftp.retrbinary(
                "RETR " + inputFile, open(tmpFile.format(inputFile), "wb").write
            )

            try:
                data = pd.read_csv(tmpFile.format(inputFile))[
                    [
                        "gene symbol",
                        "disease name",
                        "allelic requirement",
                        "gene disease pair entry date",
                        "confidence category",
                        "panel",
                    ]
                ].rename(
                    columns={
                        "gene symbol": "targetFromSourceId",
                        "disease name": "diseaseFromSource",
                        "allelic requirement": "allelicRequirements",
                        "gene disease pair entry date": "entryDate",
                        "confidence category": "confidence",
                        "panel": "studyId",
                    }
                )
                # cleaning disease name to match evidence file
                data["diseaseFromSource"] = data["diseaseFromSource"].str.replace(
                    r".+-related ", "", regex=True
                )
                # fetching year
                data["year"] = data.apply(
                    lambda row: (
                        int(row.entryDate.split("-")[0])
                        if (row.entryDate == row.entryDate)
                        else None
                    ),
                    axis=1,
                )
                data = data.drop(columns=["entryDate"])
                data["datasourceId"] = "gene2phenotype"
                data.to_csv(
                    outputFile.format(inputFile.replace(".csv.gz", ".csv")), index=False
                )

                print(
                    outputFile.format(inputFile.replace(".csv.gz", ".csv")),
                    "generated",
                )
                os.remove(tmpFile.format(inputFile))

            except:
                os.remove(outputFile.format(inputFile.replace(".csv.gz", ".csv")))
                os.remove(tmpFile.format(inputFile))
                print(inputFile, "error")

    ftp.quit()


# getYearFromPubMedBaseline()
# getYearFromPubMedUpdatefiles()
# getYearFromExpressionAtlas()
# getYearFromCancerGenomeInterpreter()
# getYearFromClinGen()
# getYearFromGene2Phenotype()

# other datasources:
#  ot_genetics_portal publicationYear is never null! :)
#  europepmc publicationYear is never null! :)
#  sysbio 100% covered by pmid -> map to PubMed years
#  eva and eva_somatic provisional solution: gs://ot-team/dsuveges/ClinVar_dated_release.2023.07
#  intogene - releases since 2013 no mention to publication date
#  orphanet - unsuccessfully tried to find dates in associations file https://www.orphadata.com/genes/
#  progeny - no date - @irene: What we call Progeny is the application of this algorithm to gene expression data from TCGA. The data is here https://storage.googleapis.com/otar000-evidence_input/PROGENy/data_files/progeny_normalVStumor_opentargets.txt
#  genomics_england - @irene's link doesn't have version year: https://storage.googleapis.com/otar000-evidence_input/GenomicsEngland/data_files/All_genes_20220804-1350_green_amber_red_public_v1plus_no_superpanels.tsv
#  cancer_gene_census
#  reactome
#  impc - no mention to date in files
# pubmed and clinvar too large to tun them locally

# Paths
ot_platform_version = "23.06"

## data path
data_path = "/Users/mariaf/OT_platform/{}/".format(ot_platform_version)
evidence_file = data_path + "evidence"
targets_file = data_path + "targets"
diseases_file = data_path + "diseases"

## results path
results_path = "/Users/mariaf/TargetEngine/results/{}/".format(ot_platform_version)
evidenceDated_file = results_path + "evidenceDated"


def addYearFromGeneBass(version="2022/04/13"):
    spark = SparkSession.builder.getOrCreate()

    data = getYearFromGeneBass(version=version)
    data = spark.createDataFrame(
        [[data[d] for d in sorted(data)]],
        schema=sorted(data),
    )

    results = (
        spark.read.parquet(evidence_file + "/sourceId=gene_burden")
        .withColumnRenamed(
            "literature", "literatureOriginal"
        )  # it's always null when projectId=Genebass
        .join(
            data.withColumnRenamed("literature", "literatureGenebass"),
            ["datasourceId", "projectId"],
            "left",
        )
        .withColumn(
            "literature",
            F.when(
                F.col("projectId") == "Genebass", F.col("literatureGenebass")
            ).otherwise(F.col("literatureOriginal")),
        )
        .drop("literatureOriginal", "literatureGenebass")
    )
    results.write.parquet(evidenceDated_file + "/sourceId=gene_burden_tmp")


def addYearFromSlapenrich():
    spark = SparkSession.builder.getOrCreate()

    data = getYearFromSlapenrich()
    data = spark.createDataFrame(
        [[data[d] for d in sorted(data)]],
        schema=sorted(data),
    )

    results = (
        spark.read.parquet(evidence_file + "/sourceId=slapenrich")
        .drop("literature")  # it's always null
        .join(
            data,
            ["datasourceId"],
            "left",
        )
    )
    results.write.parquet(evidenceDated_file + "/sourceId=slapenrich_tmp")


def addYearFromCRISPRbrain():
    spark = SparkSession.builder.getOrCreate()

    data = getYearFromCRISPRbrain()
    data = spark.createDataFrame(
        [[data[d] for d in sorted(data)]],
        schema=sorted(data),
    )

    results = (
        spark.read.parquet(evidence_file + "/sourceId=crispr_screen")
        .withColumnRenamed(
            "literature", "literatureOriginal"
        )  # sometimes it's null and sometimes it's not when projectId=crispr_brain
        .join(
            data.withColumnRenamed("literature", "literatureCrisprbrain"),
            ["datasourceId", "projectId"],
            "left",
        )
        .withColumn(
            "literature",
            F.when(
                (F.col("projectId") == "crispr_brain")
                & (F.col("literatureOriginal").isNull()),
                F.col("literatureCrisprbrain"),
            ).otherwise(F.col("literatureOriginal")),
        )
        .drop("literatureOriginal", "literatureCrisprbrain")
    )
    results.write.parquet(evidenceDated_file + "/sourceId=crispr_screen_tmp")


def addYearFromProjectScore():
    spark = SparkSession.builder.getOrCreate()

    data = getYearFromProjectScore()
    data = spark.createDataFrame(
        [[data[d] for d in sorted(data)]],
        schema=sorted(data),
    )

    results = (
        spark.read.parquet(evidence_file + "/sourceId=crispr")
        .drop("literature")  # it's always null
        .join(
            data,
            ["datasourceId"],
            "left",
        )
    )
    results.write.parquet(evidenceDated_file + "/sourceId=crispr_tmp")


def addYearFromExpressionAtlas(version=2023.06):
    # 3 experiments have forbidden access: E-MTAB-9045, E-MTAB-8887, E-MTAB-9950. They map to the undated pieces of evidence
    spark = SparkSession.builder.getOrCreate()

    expressionatlas_path = "{}/ExpressionAtlas/{}/".format(
        db_path, version
    )  # Pedro Madrigal, contact person, http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments

    results = (
        spark.read.parquet(evidence_file + "/sourceId=expression_atlas")
        .join(
            spark.read.option("header", "true")
            .csv(expressionatlas_path + "experiment2year.csv")
            .withColumnRenamed("year", "studyYear"),
            ["studyId", "datasourceId"],
            "left",
        )
        .withColumn("studyYear", F.col("studyYear").cast("long"))
    )
    results.write.parquet(evidenceDated_file + "/sourceId=expression_atlas_tmp")


# do not run this because getYearFromCancerGenomeInterpreter() function messes around with sources field and results in more evidence than originally,
# until this is fixed, just map pmids in literature to pubmedYear
def addYearFromCancerGenomeInterpreter(
    version="2018.01.17",
):  # Last update: 2018/01/17 (from website)
    spark = SparkSession.builder.getOrCreate()

    cancergenomeinterpreter_path = "{}/CancerGenomeInterpreter/{}/".format(
        db_path,
        version,
    )  # https://www.cancergenomeinterpreter.org/2018/data/cgi_biomarkers_latest.zip -> https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/cancerBiomarkers.py

    results = (
        spark.read.parquet(evidence_file + "/sourceId=cancer_biomarkers")
        .join(
            spark.read.option("header", "true")
            .csv(cancergenomeinterpreter_path + "evidence2year.csv")
            .withColumnRenamed("year", "curationYear"),
            [
                "targetFromSourceId",
                "diseaseFromSource",
                "datasourceId",
                "drugFromSource",
                "biomarkerName",
                "confidence",
            ],
            "left",
        )
        .withColumn("curationYear", F.col("curationYear").cast("long"))
    )
    results.write.parquet(evidenceDated_file + "/sourceId=cancer_biomarkers_tmp")


def addYearFromClinGen(version="2023.07.25"):
    spark = SparkSession.builder.getOrCreate()
    clingen_path = "{}/ClinGen/{}/".format(db_path, version)

    results = (
        spark.read.parquet(evidence_file + "/sourceId=clingen")
        .withColumn(
            "url", F.explode_outer("urls.url")
        )  # explode results in only one url per urls
        .join(
            spark.read.option("header", "true")
            .csv(clingen_path + "evidence2year.csv")
            .withColumnRenamed("year", "curationYear"),
            [
                "targetFromSourceId",
                "diseaseFromSourceId",
                "url",  # never null
                "confidence",
                "studyId",
                "datasourceId",
            ],
            "left",
        )
        .drop("url")
        .withColumn("curationYear", F.col("curationYear").cast("long"))
    )
    results.write.parquet(evidenceDated_file + "/sourceId=clingen_tmp")


def addYearFromGene2Phenotype(version="2023.07.25"):
    # some pairs lack entry date in raw data
    spark = SparkSession.builder.getOrCreate()
    gene2phenotype_path = "{}/Gene2Phenotype/{}/".format(
        db_path, version
    )  # http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/28_04_2023/

    results = (
        spark.read.parquet(evidence_file + "/sourceId=gene2phenotype")
        .withColumn(
            "tmp", F.explode_outer("allelicRequirements")
        )  # explode results in only one allelicRequirement per allelicRequirements
        .fillna("0", subset=["tmp"])
        .join(
            spark.read.option("header", "true")
            .csv(gene2phenotype_path + "pair2year")
            .withColumnRenamed("allelicRequirements", "tmp")
            .fillna("0", subset=["tmp"])
            .withColumnRenamed("year", "curationYear"),
            [
                "targetFromSourceId",
                "diseaseFromSource",
                "tmp",
                "confidence",
                "studyId",  # is never null
                "datasourceId",
            ],
            "left",
        )
        .drop("tmp")
        .withColumn("curationYear", F.col("curationYear").cast("long"))
    )
    results.write.parquet(evidenceDated_file + "/sourceId=gene2phenotype_tmp")


addYearFromGeneBass()
addYearFromSlapenrich()
addYearFromCRISPRbrain()
addYearFromProjectScore()
addYearFromExpressionAtlas()
addYearFromClinGen()
addYearFromGene2Phenotype()
