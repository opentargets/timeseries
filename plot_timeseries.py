#!/usr/bin/env python

__author__ = "Maria J. Falaguera (mjfalagueramata@gmail.com)"
__date__ = "02 Jun 2025"

"""
plot_timeries.py: Plot timeseries for associations in Open Targets Platform.
"""

import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
import seaborn as sns
import sys
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# Settings
font = {"weight": "normal", "size": 16, "family": "Arial"}
palette = "tab20"
overallColor = "darkgray"
geneticColor = "crimson"
literatureColor = "dodgerblue"
clinicalColor = "limegreen"
firstYear = 1995
lastYear = 2023
ymin = 0
ymax = 1
transparent = False
figsizeLandscape = (6, 3.5)
figsizeSquare = (3.5, 3.5)
cutoffLinewidth = 0.1
thickLineWidth = 5
slimLineWidth = 1
alpha = 0.8
evidenceAlpha = 0.2
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
    {
        "id": "gwas_credible_sets",
        "sectionId": "gwas_credible_sets",
        "label": "GWAS associations",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gwas-associations",
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
datatypeOrder = [
    "Known drug",
    "Animal model",
    "Literature",
    "Genetic association",
    "Somatic mutations",
    "RNA expression",
    "Affected pathway",
    "All",
]

datasourceOrder = [
    "Gene2phenotype",
    "Cancer Biomarkers",
    "ClinVar (somatic)",
    "Clingen",
    "ClinVar",
    "IMPC",
    "Reactome",
    "UniProt curated variants",
    "UniProt literature",
    "Orphanet",
    "ChEMBL",
    "GEL PanelApp",
    "Expression Atlas",
    "Cancer Gene Census",
    "Gene signatures",
    "SLAPenrich",
    "Europe PMC",
    "CRISPR Screens",
    "Gene Burden",
    "GWAS associations",
    "OT Genetics",
    "Project Score",
]


sns.set(font_scale=1.5)
sns.set_style("ticks")

# Paths
ot_platform_version = "/23.06/"

## data path
data_path = "/Users/mariaf/OT_platform/{}/".format(ot_platform_version)
targets_file = data_path + "targets"
diseases_file = data_path + "diseases"

## results path
results_path = "/Users/mariaf/TargetEngine/results/{}/".format(ot_platform_version)
evidenceIndirectDated_file = results_path + "/evidenceIndirectDated"
evidenceDirectDated_file = results_path + "/evidenceDirectDated"
associationByDatasourceIndirectOverYears_file = (
    results_path + "/associationByDatasourceIndirectOverYears"
)
associationByOverallIndirectOverYears_file = (
    results_path + "/associationByOverallIndirectOverYears"
)
associationByDatasourceDirectOverYears_file = (
    results_path + "/associationByDatasourceDirectOverYears"
)
associationByOverallDirectOverYears_file = (
    results_path + "associationByOverallDirectOverYears"
)
targetByDatasourceIndirectOverYears_file = (
    results_path + "targetByDatasourceIndirectOverYears"
)
targetByOverallIndirectOverYears_file = (
    results_path + "targetByOverallIndirectOverYears"
)

## plots path
plots_path = results_path + "plots/"


def getDatasourceToName(output="list"):
    """
    Returns list of data sources names/labels.
    """

    id2label_dict = {
        datasource["id"]: datasource["label"] for datasource in dataSources
    }
    id2label_dict["other"] = "Other sources"
    if output == "dict":
        return id2label_dict
    elif output == "list":
        id2label_list = [[idx, id2label_dict[idx]] for idx in id2label_dict]
        return id2label_list


def getDatatypeForDatasource(output="list"):
    """
    Returns list of data sources names/labels.
    """

    source2type_dict = {
        datasource["id"]: datasource["aggregationId"] for datasource in dataSources
    }
    if output == "dict":
        return source2type_dict
    elif output == "list":
        source2type_list = [
            [source, source2type_dict[source]] for source in source2type_dict
        ]
        return source2type_list
    elif output == "dataframe":
        return pd.DataFrame(dataSources).rename(
            columns={"aggregationId": "datatypeId", "id": "datasourceId"}
        )[["datatypeId", "datasourceId"]]


def getDatasourceToColour1(by="label", output="dict", palette="tab20", order=False):
    """
    Returns colourmap dictionary of data sources.
    """

    colours = {}

    if by == "label":
        if order:
            labels = order
        else:
            labels = sorted(list(getDatasourceToName(output="dict").values()))

        for label, colour in zip(
            labels,
            sns.color_palette(palette=palette, n_colors=len(labels)).as_hex(),
        ):
            colours[label] = colour

        colours["Cancer Gene Census"] = "#9edae5"  # light blue
        colours["Cancer Biomarkers"] = "#17becf"  # turquoise
        colours["ClinVar"] = "#e377c2"  # pink
        colours["ClinVar (somatic)"] = "#f7b6d2"  # light pink
        colours["Europe PMC"] = "#bcbd22"  # light green
        colours["Orphanet"] = "#dbdb8d"  # very light green
        colours["ChEMBL"] = "#d62728"  # red
        colours["Clingen"] = "#ff9896"  # light red
        colours["UniProt literature"] = "#c5b0d5"  # light purple
        colours["UniProt curated variants"] = "#9467bd"  # purple
        colours["IntOGen"] = "#2ca02c"  # green
        colours["SLAPenrich"] = "#98df8a"  # light green
        colours["All"] = "#D3D3D3"
        colours["OT Genetics"] = "#1f77b4"  # blue
        colours["IMPC"] = "#ff7f0e"  # orange
        colours["Gene Burden"] = "#8c564b"  # brown
        colours["Gene2phenotype"] = "#c49c94"  # light brown
        colours["GEL PanelApp"] = "#2ca02c"  # dark green
        colours["Project Score"] = "#ffbb78"  # light orange

    elif by == "id":
        ids = sorted(list(getDatasourceToName(output="dict").keys()))

        for idx, colour in zip(
            ids,
            sns.color_palette(palette=palette, n_colors=len(ids)).as_hex(),
        ):
            colours[idx] = colour

        colours["cancer_gene_census"] = "#9edae5"  # light blue
        colours["cancer_biomarkers"] = "#17becf"  # turquoise
        colours["eva"] = "#e377c2"  # pink
        colours["eva_somatic"] = "#f7b6d2"  # light pink
        colours["europepmc"] = "#bcbd22"  # light green
        colours["orphanet"] = "#dbdb8d"  # very light green
        colours["chembl"] = "#d62728"  # red
        colours["clingen"] = "#ff9896"  # light red
        colours["uniprot_literature"] = "#c5b0d5"  # light purple
        colours["uniprot_curated_variants"] = "#9467bd"  # purple
        colours["intogen"] = "#2ca02c"  # green
        colours["slapenrich"] = "#98df8a"  # light green
        colours["all"] = "#D3D3D3"
        colours["ot_genetics_portal"] = "#1f77b4"  # blue
        colours["gwas_credible_sets"] = "#1f77b4"  # blue
        colours["impc"] = "#ff7f0e"  # orange
        colours["gene_burden"] = "#8c564b"  # brown
        colours["gene2phenotype"] = "#c49c94"  # light brown
        colours["genomics_england"] = "#2ca02c"  # dark green
        colours["crispr"] = "#ffbb78"  # light orange

    if output == "dict":
        return colours
    elif output == "list":
        return [[source, colour] for source, colour in colours.items()]


def getDatasourceToColour():
    result = {
        datasourceId: getDatatypeToColour()[datatypeId]
        for datasourceId, datatypeId in getDatatypeForDatasource(output="dict").items()
    }
    result["other"] = "lightgray"
    return result


def getDatasourceToShape(by="label", output="dict", order=False):
    """
    Returns shape dictionary of data sources.
    """

    shapes = {}

    if by == "label":
        if order:
            labels = order
        else:
            labels = sorted(list(getDatasourceToName(output="dict").values()))

        for label in labels:
            if label == "Europe PMC":
                shapes[label] = "o"
            elif label == "ChEMBL":
                shapes[label] = "o"
            else:
                shapes[label] = "o"

    elif by == "id":
        ids = sorted(list(getDatasourceToName(output="dict").keys()))

        for idx in ids:
            if idx == "europepmc":
                shapes[label] = "o"
            elif label == "chembl":
                shapes[label] = "o"
            else:
                shapes[label] = "o"

    if output == "dict":
        return shapes
    elif output == "list":
        return [[source, shape] for source, shape in shapes.items()]


def getDatatypeToName(output="list"):
    id2label_dict = {
        datasource["aggregationId"]: datasource["aggregation"]
        for datasource in dataSources
    } | {
        "clinical": "Clinical",
        "somatic_genetic": "Human genetic data",
        "0.5": "Pre-clinical phase",
        "1.0/2.0": "Clinical phase I/II",
        "3.0": "Clinical phase III",
        "4.0": "Clinical phase IV",
        "other": "Other",
        "europepmc": "Literature",
    }

    if output == "dict":
        return id2label_dict
    elif output == "list":
        id2label_list = [[idx, id2label_dict[idx]] for idx in id2label_dict]
        return id2label_list


def getDatatypeToColour(
    by="id",
):
    """
    Returns colourmap dictionary of data sources.
    """

    colours = {
        "genetic_association": "#4280BE",
        "clinical": "#FF595E",
        "known_drug": "#FF595E",
        "literature": "#AEB236",
        "animal_model": "#fed154",  # "#E8C468",
        "somatic_mutation": "#9CC5E1",
        "rna_expression": "#63CCCA",
        "affected_pathway": "#f4a261",
        "somatic_genetic": "#4280BE",
        "4.0": "#F2C8BE",
        "3.0": "#FC8F92",
        "1.0/2.0": "#FF595E",  # "#e76f51",
        "0.5": "gray",
        "other": "lightgray",
    }

    if by == "label":
        return {
            getDatatypeToName(output="dict")[datatypeId]: colours[datatypeId]
            for datatypeId in colours
        }

    else:
        return colours


def getTherapeuticAreaToName():
    id2name = {
        "EFO_0000319": "cardiovascular disease",
        "EFO_0000540": "immune system disease",
        "EFO_0000618": "nervous system disease",
        "EFO_0000651": "phenotype",
        "EFO_0001379": "endocrine system disease",
        "EFO_0005741": "infectious disease",
        "EFO_0005803": "hematologic disease",
        "EFO_0009605": "pancreas disease",
        "EFO_0009690": "urinary system disease",
        "EFO_0010282": "gastrointestinal disease",
        "EFO_0010285": "integumentary system disease",
        "MONDO_0002025": "psychiatric disorder",
        "MONDO_0021205": "disorder of ear",
        "MONDO_0024458": "disorder of visual system",
        "MONDO_0045024": "cancer or benign tumor",
        "OTAR_0000006": "musculoskeletal or connective tissue disease",
        "OTAR_0000009": "injury, poisoning or other complication",
        "OTAR_0000010": "respiratory or thoracic disease",
        "OTAR_0000014": "pregnancy or perinatal disease",
        "OTAR_0000017": "reproductive system or breast disease",
        "OTAR_0000018": "genetic, familial or congenital disease",
        "OTAR_0000020": "nutritional or metabolic disease",
    }
    return id2name


def getTherapeuticAreaToColour(
    by="name",
    output="dict",
    palette="tab10",
):
    """
    Returns colourmap dictionary of data sources.
    """

    colours = {}

    if by == "name":
        names = list(getTherapeuticAreaToName().values())

        for name, colour in zip(
            names,
            sns.color_palette(palette=palette, n_colors=len(names)).as_hex()[::-1],
        ):
            colours[name] = colour

    elif by == "id":
        ids = list(getTherapeuticAreaToName().keys())

        for idx, colour in zip(
            ids,
            sns.color_palette(palette=palette, n_colors=len(ids)).as_hex()[::-1],
        ):
            colours[idx] = colour

    colours["All"] = "#D3D3D3"

    if output == "dict":
        return colours
    elif output == "list":
        return [[name, colour] for name, colour in colours.items()]


def plotTargetDisease(
    targetId,
    diseaseId,
    includeNonDated=True,
    showLegend=True,
    showScore=False,
    showEvidence=True,
    showNovelty=True,
    showOverall=True,
    img=None,
    transparent=transparent,
    vlines=[],
    dashedvlines=[],
    hlines=[],
    annotate=False,
    evidenceLink="indirect",
    title=True,
    xlabel=False,
    ylabel="Score",
    area=False,
    xlim=None,
    figsize="landscape",
):
    """
    Plot timeseries for a disease-target associations.

    Args:
        targetId (str):         target
        diseaseId (str):        disease
        includeNonDated (bool)  include non dated evidence labeled as lastYear + 1
        major (str):            score/novelty predominant curve
        minor (str):            score/novelty secondary curve
        img (str):              path to save image
        vlines (list):          verical lines to plot
        hlines ( list):         horizontal lines to plot

    Returns:
        img with timeseries plot

    """

    # Establish spark connection
    spark = SparkSession.builder.getOrCreate()

    if evidenceLink == "indirect":
        overall_data = associationByOverallIndirectOverYears_file
        datasource_data = associationByDatasourceIndirectOverYears_file
        evidence_data = evidenceIndirectDated_file
    elif evidenceLink == "direct":
        overall_data = associationByOverallDirectOverYears_file
        datasource_data = associationByDatasourceDirectOverYears_file
        evidence_data = evidenceDirectDated_file

    # overall data
    overall_data = (
        spark.read.parquet(overall_data)
        # filter target and disease
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
        # add target name
        .join(
            spark.read.parquet(targets_file).select(
                F.col("id").alias("targetId"),
                F.col("approvedSymbol").alias("targetSymbol"),
            ),
            "targetId",
            "left",
        )
        # add disease name
        .join(
            spark.read.parquet(diseases_file).select(
                F.col("id").alias("diseaseId"),
                F.col("name").alias("diseaseName"),
            ),
            "diseaseId",
            "left",
        ).orderBy("year")
    )

    # datasource data
    datasourceNames = spark.createDataFrame(
        data=[
            [datasourceId, datasourceName]
            for datasourceId, datasourceName in getDatasourceToName()
        ],
        schema=["datasourceId", "datasourceName"],
    )
    datasource_data = (
        spark.read.parquet(datasource_data)
        # filter disease and target
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
        # add datasources' names
        .join(datasourceNames, "datasourceId", "left")
        .orderBy("datasourceName")
        .orderBy("year")
    )
    # evidence data
    evidence_data = (
        spark.read.parquet(evidence_data)
        # filter disease and target
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
        # add datasources' names
        .join(datasourceNames, "datasourceId", "left")
        .orderBy("datasourceName")
        .orderBy("year")
    )

    if not includeNonDated:
        overall_data = overall_data.join(
            overall_data.filter(F.col("year").isNotNull() & (F.col("score") > 0))
            .select("diseaseId", "targetId")
            .distinct(),
            ["diseaseId", "targetId"],
            "inner",
        )
        datasource_data = datasource_data.join(
            datasource_data.filter(F.col("year").isNotNull() & (F.col("score") > 0))
            .select("datasourceId")
            .distinct(),
            ["datasourceId"],
            "inner",
        )
        evidence_data = evidence_data.join(
            evidence_data.filter(F.col("year").isNotNull() & (F.col("score") > 0))
            .select("datasourceId")
            .distinct(),
            ["datasourceId"],
            "inner",
        )

    # pandas
    overall_data = overall_data.toPandas()
    datasource_data = datasource_data.toPandas()
    evidence_data = evidence_data.toPandas()

    # include non-dated evidence labeled as lastYear + 1
    if includeNonDated:
        overall_data.year = overall_data.year.fillna(lastYear + 1)
        datasource_data.year = datasource_data.year.fillna(lastYear + 1)
        evidence_data.year = evidence_data.year.fillna(lastYear + 1)

    # initialize figure
    plt.figure()
    if figsize == "landscape":
        figsize = figsizeLandscape
    elif figsize == "square":
        figsize = figsizeSquare
    fig, ax = plt.subplots(figsize=figsize)
    matplotlib.rc("font", **font)
    plt.rcParams["savefig.facecolor"] = "white"

    # extra features on the background
    for vline in vlines:
        print(vline)
        ax.axvline(x=vline, lw=1, linestyle="-", color="k")
    for vline in dashedvlines:
        print(vline)
        ax.axvline(x=vline, lw=1, linestyle="--", color="k")

    for hline in hlines:
        ax.axhline(y=hline, lw=1, linestyle="-", color="k")
    if area:
        ax.axvspan(area[0], area[1], alpha=0.1, color="k")

    # blank profile for cases with no evidence
    sns.scatterplot(
        data=pd.DataFrame(
            [
                [
                    2010,
                    0.0,
                    "a",
                ]
            ],
            columns=["year", "score", "datasourceName"],
        ),
        x="year",
        y="score",
        hue="datasourceName",
        style="datasourceName",
        s=200,
        alpha=0,
        legend=False,
        ax=ax,
    )

    if showEvidence:

        # scatter plot of evidence
        sns.scatterplot(
            data=evidence_data,
            x="year",
            y="score",
            hue="datasourceName",
            style="datasourceName",
            palette=dict(getDatasourceToColour()),
            markers=dict(getDatasourceToShape()),
            s=200,
            alpha=evidenceAlpha,
            ax=ax,
        )

    if showScore:

        # line plot of source score
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y="score",
            hue="datasourceName",
            palette={source: "w" for source in dict(getDatasourceToColour())},
            lw=slimLineWidth + 4,
            markersize=0,
            linestyle="-",
            ax=ax,
            alpha=1,
        )
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y="score",
            hue="datasourceName",
            palette=dict(getDatasourceToColour()),
            lw=slimLineWidth,
            markersize=0,
            linestyle="-",
            ax=ax,
            alpha=0.8,
        )

        if showOverall:

            # line plot of overall score
            sns.lineplot(
                data=overall_data.fillna(0),
                x="year",
                y="score",
                color="w",
                lw=slimLineWidth + 4,
                marker="o",
                markersize=0,
                linestyle="-",
                ax=ax,
                alpha=1,
            )
            sns.lineplot(
                data=overall_data.fillna(0),
                x="year",
                y="score",
                color="k",
                lw=slimLineWidth,
                marker="o",
                markersize=0,
                linestyle="--",
                ax=ax,
                alpha=0.8,
            )

    if showNovelty:

        # line plot of source novelty
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y="novelty",
            hue="datasourceName",
            palette={source: "w" for source in dict(getDatasourceToColour())},
            lw=thickLineWidth + 4,
            markersize=0,
            linestyle="-",
            ax=ax,
            legend=True,
            alpha=1,
        )
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y="novelty",
            hue="datasourceName",
            palette=dict(getDatasourceToColour()),
            lw=thickLineWidth,
            markersize=0,
            linestyle="-",
            ax=ax,
            alpha=0.8,
        )

        if showOverall:
            # line plot of overall novelty
            sns.lineplot(
                data=overall_data.fillna(0),
                x="year",
                y="novelty",
                color="w",
                lw=thickLineWidth + 4,
                marker="o",
                markersize=0,
                linestyle="-",
                ax=ax,
                alpha=1,
            )
            sns.lineplot(
                data=overall_data.fillna(0),
                x="year",
                y="novelty",
                color="k",
                lw=thickLineWidth - 1,
                marker="o",
                markersize=0,
                linestyle="--",
                ax=ax,
                alpha=0.8,
            )

    if showLegend:
        legend_elements = []
        for source in sorted(datasource_data.datasourceName.unique()):
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker=getDatasourceToShape()[source],
                    color="w",
                    label=source + " evidence",
                    markerfacecolor=dict(getDatasourceToColour())[source],
                    alpha=evidenceAlpha + 0.1,
                    markersize=20,
                )
            )
        for source in sorted(datasource_data.datasourceName.unique()):
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=dict(getDatasourceToColour())[source],
                    lw=slimLineWidth,
                    label=source + " score",
                )
            )
        if showOverall:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=overallColor,
                    lw=slimLineWidth,
                    label="Overall score",
                )
            )
        for source in sorted(datasource_data.datasourceName.unique()):
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=dict(getDatasourceToColour())[source],
                    lw=thickLineWidth,
                    label=source + " novelty",
                )
            )
        if showOverall:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=overallColor,
                    lw=thickLineWidth,
                    label="Overall novelty",
                )
            )

        ax.legend(
            handles=legend_elements,
            frameon=False,
            fontsize=font["size"],
            loc="center left",
            bbox_to_anchor=(1, 0.5),
        )
    else:
        try:
            ax.get_legend().remove()
        except AttributeError:
            print("No data")

    # title
    if title == True:
        ax.set_title(
            "{} + {}".format(
                overall_data.diseaseName.values[0].capitalize(),
                overall_data.targetSymbol.values[0],
            ),
            fontfamily=font["family"],
            fontsize=font["size"],
        )
    elif title == False:
        ax.set_title("")
    else:
        ax.set_title(title)

    # xlabel
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel("")

    # customize plot
    ax.set_ylabel(
        ylabel,
        # fontsize=font["size"],
    )
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    elif includeNonDated:
        ax.set_xlim(firstYear, lastYear + 1)
        ax.text(lastYear - 12, 0.9, "Not dated:")
    else:
        ax.set_xlim(firstYear, lastYear)
    ax.set_ylim(ymin, ymax)

    # extra features
    if includeNonDated:
        ax.axvline(x=lastYear, lw=1, linestyle="-", color="k")
    if annotate:
        ax.annotate(annotate, (1995.5, 0.9))
        # ax.annotate(annotate, (2000.5, 0.9))
    if img is not None:
        fig.savefig(img, bbox_inches="tight", transparent=transparent, dpi=300)
        plt.close()
    ax.axhline(y=0, lw=1, linestyle="-", color="k", alpha=0)

    overall_data["datasourceId"] = "overall"
    overall_data["sourceId"] = "overall"
    overall_data["datasourceName"] = "Overall"
    return pd.concat([datasource_data.fillna(0), overall_data.fillna(0)])
