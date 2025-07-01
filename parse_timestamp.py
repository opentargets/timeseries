#!/usr/bin/env python

__author__ = "Cote Falaguera (mjfalagueramata@gmail.com)"
__date__ = "21 May 2025"

import xml.etree.ElementTree as ET
import os
import csv
from datetime import datetime, timedelta
import ftplib
import pandas as pd
import requests
import inspect
import subprocess
import time
from time import sleep
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
from concurrent.futures import ThreadPoolExecutor, as_completed


"""
parse_timestamp.py: Parse timestamp data from raw sources for Open Targets Platform evidence and upload them to GCS.

Useful GitHub links:
- https://github.com/opentargets/timeseries
- https://github.com/opentargets/issues/issues/2739
"""

# Paths
timestamp_version = datetime.now().strftime("%d.%m.%Y")
ot_version = "25.03"
timestamp_path = "/Users/mariaf/TargetEngine/data/timestamp/{}/".format(
    timestamp_version
)
evidence_path = "/Users/mariaf/OT_platform/{}/evidence".format(ot_version)
gcs_path = "gs://ot-team/cfalaguera/"


# Functions to retrieve timestamps for undated evidence from original sources


## Expression Atlas
def fetch_release_year_from_idf(
    accession, max_retries=3
):  # accession corresponds to studyId field in evidence/source=expression_atlas file

    # FTP URL pattern
    # I query individual experiments to avoid downloading the whole database which takes longer
    url = f"https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/{accession}/{accession}.idf.txt"
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=20)
            response.raise_for_status()
            for line in response.text.splitlines():
                if line.startswith("Public Release Date"):
                    parts = line.strip().split("\t")
                    if len(parts) > 1:
                        date = parts[1]
                        year = date[:4]
                        return accession, year
            # If not found, return empty
            return accession, ""
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"{accession}: error ({e}), retrying...")
            else:
                print(f"{accession}: failed after {max_retries} attempts ({e})")
    return accession, ""


def fetch_all_release_years(accessions, output, max_workers=16):
    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["studyId", "year"])

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_accession = {
                executor.submit(fetch_release_year_from_idf, acc): acc
                for acc in accessions
            }
            total = len(accessions)
            processed = 0
            for future in as_completed(future_to_accession):
                accession, release_year = future.result()
                writer.writerow([accession, release_year])
                processed += 1
                if processed % 1000 == 0 or processed == total:
                    print(
                        f"Processed {processed}/{total} accessions ({processed/total:.1%})"
                    )


def get_year_from_expression_atlas():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    expression_atlas_path = timestamp_path + "/expression_atlas/"
    os.makedirs(expression_atlas_path + "/tmp/", exist_ok=True)

    with SparkSession.builder.getOrCreate() as spark:
        fetch_all_release_years(
            accessions=spark.read.parquet(evidence_path + "/sourceId=expression_atlas")
            .select("studyId")
            .distinct()
            .toPandas()
            .studyId.tolist(),
            output=expression_atlas_path + "/evidence2year.csv",
            max_workers=16,
        )

    os.rmdir(expression_atlas_path + "/tmp/")


## Gene2phenotype
def get_year_from_gene2phenotype():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    # Ideally, I'd like to use the field gene disease pair entry date available in http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/ and extracted in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/Gene2Phenotype.py
    # but this field is not provided in the evidence file yet. In the meantime, I fetch the field myself from http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/
    # note that some of the evidence have pmid annotated and some other don't but all of them have "gene disease pair entry date", thus, I choose to use the date field instead of the pmid for the timestamp.
    url = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd("pub/databases/gene2phenotype/G2P_data_downloads/")

    # Find latest release directory
    folders = []
    ftp.retrlines("LIST", folders.append)
    last_release = max([line.split()[-1] for line in folders if line.startswith("d")])
    ftp.cwd(last_release)

    gene2phenotype_path = timestamp_path + "/gene2phenotype/"
    os.makedirs(gene2phenotype_path + "/tmp/", exist_ok=True)

    header = True
    for inputFile in ftp.nlst("*.csv.gz"):
        ftp.retrbinary(
            "RETR " + inputFile,
            open(gene2phenotype_path + "/tmp/" + inputFile, "wb").write,
        )

        try:
            data = pd.read_csv(gene2phenotype_path + "/tmp/" + inputFile)[
                [
                    "gene symbol",
                    "disease name",
                    "allelic requirement",
                    "date of last review",
                    "confidence",
                    "panel",
                ]
            ].rename(
                columns={
                    "gene symbol": "targetFromSourceId",
                    "disease name": "diseaseFromSource",
                    "allelic requirement": "allelicRequirements",
                    "date of last review": "entryDate",
                    "confidence": "confidence",
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
                gene2phenotype_path + "/evidence2year.csv",
                mode="a",
                index=False,
                header=header,
            )
            os.remove(gene2phenotype_path + "/tmp/" + inputFile)

        except:
            os.remove(gene2phenotype_path + "/tmp/" + inputFile.replace(".gz", ""))
            print(inputFile, "error")

        header = False

    os.rmdir(gene2phenotype_path + "/tmp/")
    ftp.quit()


## Cancer Genome Interpreter
def get_year_from_cancer_genome_interpreter():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    # ideally, add the field "Curation date" to the list of fields that are already fetched from https://www.cancergenomeinterpreter.org/2018/data/cgi_biomarkers_latest.zip in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/cancerBiomarkers.py
    cancer_genome_interpreter_path = timestamp_path + "/cancer_genome_interpreter/"
    os.makedirs(cancer_genome_interpreter_path + "/tmp/", exist_ok=True)

    # in the meantime, I'm gonna map the evidence in the OT file to their Curation date in the cancer biomarkers source file. For that, I need to use the drug info as a reference.
    url = "https://www.cancergenomeinterpreter.org/2018/data/cgi_biomarkers_latest.zip"
    response = requests.get(url)

    # Set the destination file path
    zip_file = os.path.join(cancer_genome_interpreter_path, "cgi_biomarkers_latest.zip")

    # Download and save the ZIP file
    url = "https://www.cancergenomeinterpreter.org/2018/data/cgi_biomarkers_latest.zip"
    response = requests.get(url)
    if response.status_code == 200:
        with open(zip_file, "wb") as f:
            f.write(response.content)
    else:
        print("Failed to download {}".format(zip_file))

    data = pd.read_csv(zip_file, sep="\t")[
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
    )  # in this case, only semicolon-separated strings have been found, which may contain inside commas that we don't want to remove
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
    # hence, I ignore this field
    data = data.rename(
        columns={
            "Gene": "targetFromSourceId",
            "Primary Tumor type full name": "diseaseFromSource",
            "Curation date": "curationDate",
            "Biomarker": "biomarkerName",
            "Evidence level": "confidence",
        }
    ).drop(columns=["Drug full name", "Drug", "Source"])
    # fetch year from curationDate field (which can take different formats: dot-separated date, slash-separated date, only last 2 digist of the year, ...)
    data["year"] = data.apply(
        lambda row: (
            None
            if pd.isna(row["curationDate"])
            else (
                int(row["curationDate"].split(".")[-1])
                if "." in row["curationDate"]
                else (
                    int("20" + row["curationDate"].split("/")[-1])
                    if "/" in row["curationDate"]
                    else None
                )
            )
        ),
        axis=1,
    )

    data = data.drop(columns=["curationDate"])
    data["datasourceId"] = "cancer_biomarkers"
    data.to_csv(cancer_genome_interpreter_path + "evidence2year.csv", index=False)
    os.remove(zip_file)
    os.rmdir(cancer_genome_interpreter_path + "/tmp/")


## ClinGen
def get_year_from_clingen():
    print("running {} ...".format(inspect.currentframe().f_code.co_name))
    # ideally, use CLASSIFICATION DATE field in https://search.clinicalgenome.org/kb/gene-validity/download which is already extracted in https://github.com/opentargets/evidence_datasource_parsers/blob/master/modules/ClinGen.py
    # in the meantime, I download https://search.clinicalgenome.org/kb/gene-validity/download and fetch the field myself.
    # watch out! 27 very recently curated associations are missing in the raw file but appear in the website I don't know why

    clingen_path = timestamp_path + "/clingen/"
    os.makedirs(clingen_path + "/tmp/", exist_ok=True)

    # Set the target file path (the file name is 'download' as per the URL)
    file_path = os.path.join(clingen_path, "Clingen-Gene-Disease-Summary.tsv")

    # Download the file
    url = "https://search.clinicalgenome.org/kb/gene-validity/download"
    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "wb") as f:
            f.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

    # Find the line number where the real header starts
    with open(file_path, encoding="utf-8") as f:
        for i, line in enumerate(f):
            if line.startswith('"GENE SYMBOL"'):
                header_line = i
                break

    # Now read the file, skipping all lines before the header
    data = pd.read_csv(
        file_path,
        sep=",",
        header=header_line,
        skip_blank_lines=True,
        dtype=str,  # optional: treat all as string
        quotechar='"',
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
    data = data[~data["year"].str.contains("\+\+", na=False)]
    data["year"] = data.apply(lambda row: row.year.split("-")[0], axis=1)
    data["datasourceId"] = "clingen"
    data.to_csv(clingen_path + "evidence2year.csv", index=False)
    os.remove(file_path)
    os.rmdir(clingen_path + "/tmp/")


## PubMed (query only PMIDs of evidence in the Platform cause downloading the whole PubMed metada takes ages especially after the merge (v25.03))
def post_pmids(pmids, api_key=None):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi"
    data = {"db": "pubmed", "id": ",".join(pmids)}
    if api_key:
        data["api_key"] = api_key
    response = requests.post(url, data=data)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    webenv = root.findtext("WebEnv")
    query_key = root.findtext("QueryKey")
    return webenv, query_key


def fetch_pubmed_batch(webenv, query_key, retstart=0, retmax=1000, api_key=None):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "query_key": query_key,
        "WebEnv": webenv,
        "retstart": retstart,
        "retmax": retmax,
        "retmode": "xml",
    }
    if api_key:
        params["api_key"] = api_key
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.content


def parse_pubmed_xml(xml_content):
    results = {}
    root = ET.fromstring(xml_content)
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID")
        year = None
        year_elem = article.find(".//PubDate/Year")
        if year_elem is not None:
            year = year_elem.text
        else:
            medline_date_elem = article.find(".//PubDate/MedlineDate")
            if medline_date_elem is not None:
                year = medline_date_elem.text[:4]
        results[pmid] = year
    return results


def chunked(iterable, size):
    """Yield successive chunks from iterable of length size."""
    for i in range(0, len(iterable), size):
        yield iterable[i : i + size]


def fetch_pubmed_years(pmids, output, api_key=None):
    BATCH_SIZE = 1000
    EPOST_CHUNK_SIZE = 10000  # NCBI hard limit
    RATE_LIMIT = 0.34 if not api_key else 0.1

    all_results = {}

    print(f"Total PMIDs: {len(pmids)}")
    for chunk_idx, pmid_chunk in enumerate(chunked(pmids, EPOST_CHUNK_SIZE)):
        print(
            f"\nUploading chunk {chunk_idx+1} ({len(pmid_chunk)} PMIDs) to NCBI EPost..."
        )
        webenv, query_key = post_pmids(pmid_chunk, api_key)
        print(f"Received WebEnv and QueryKey for chunk {chunk_idx+1}.")

        for retstart in range(0, len(pmid_chunk), BATCH_SIZE):
            print(
                f"  Fetching records {retstart + 1} to {min(retstart + BATCH_SIZE, len(pmid_chunk))} in chunk {chunk_idx+1}..."
            )
            xml_content = fetch_pubmed_batch(
                webenv, query_key, retstart=retstart, retmax=BATCH_SIZE, api_key=api_key
            )
            batch_results = parse_pubmed_xml(xml_content)
            all_results.update(batch_results)
            sleep(RATE_LIMIT)

    # Write to CSV
    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pmid", "year"])
        for pmid in pmids:
            writer.writerow([pmid, all_results.get(pmid, "")])


def get_year_from_pubmed():
    pubmed_path = timestamp_path + "/pubmed/"
    os.makedirs(pubmed_path + "/tmp/", exist_ok=True)
    with SparkSession.builder.getOrCreate() as spark:

        fetch_pubmed_years(
            spark.read.parquet(evidence_path)
            .filter(F.col("publicationYear").isNull() & F.col("literature").isNotNull())
            .select(
                F.explode("literature").alias("pmid")
            )  # PubMed is huge, query only pmids with no publication year annotated yet in the evidence file "publicationYear" field
            .distinct()
            .toPandas()
            .pmid.tolist(),
            output=pubmed_path + "/evidence2year.csv",
        )

    os.rmdir(pubmed_path + "/tmp/")


# Run retriever functions
start_time = time.perf_counter()

get_year_from_gene2phenotype()
get_year_from_cancer_genome_interpreter()
get_year_from_clingen()
get_year_from_expression_atlas()
get_year_from_pubmed()


end_time = time.perf_counter()
elapsed_seconds = end_time - start_time
elapsed_td = timedelta(seconds=elapsed_seconds)
days = elapsed_td.days
hours, remainder = divmod(elapsed_td.seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"\nElapsed time: {days:02d}-{hours:02d}:{minutes:02d}:{seconds:02d}")

# Upload output files to GCS
command = [
    "gsutil",
    "cp",
    "-r",
    timestamp_path.replace("/{}/".format(timestamp_version), ""),
    gcs_path,
]
result = subprocess.run(command, capture_output=True, text=True)
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)
print("Return code:", result.returncode)
print(timestamp_path, " successfully uploaded to ", gcs_path)
