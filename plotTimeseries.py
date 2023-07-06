#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "6 Jul 2023"

"""
plotTimeries.py: Plot timeseries for associations in Open Targets Platform.
"""

import matplotlib
from matplotlib import pyplot as plt
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
import seaborn as sns


# Settings
font = {"weight": "normal", "size": 15, "family": "Verdana"}
overallColor = "lightgray"
geneticColor = "crimson"
literatureColor = "dodgerblue"
clinicalColor = "limegreen"
firstYear = 1970
lastYear = 2023
ymin = 0
ymax = 1
figsize = (15, 4)
cutoffLinewidth = 0.1
thickLineWidth = 2.5
slimLineWidth = 0.5
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
ot_platform_version = "/23.06/"

## data path
data_path = "/Users/mariaf/TargetEngine/data/OT_platform/{}/".format(
    ot_platform_version
)
targets_file = data_path + "targets"
diseases_file = data_path + "diseases"

## results path
results_path = "/Users/mariaf/TargetEngine/results/{}/".format(ot_platform_version)
associationByDatasourceIndirectOverYears_file = (
    results_path + "associationByDatasourceIndirectOverYears"
)
associationByOverallIndirectOverYears_file = (
    results_path + "associationByOverallIndirectOverYears"
)

associationByDatasourceIndirectOverYearsSignature_file = (
    results_path + "associationByDatasourceIndirectOverYearsSignature"
)
associationByLiteratureIndirectOverYears_file = (
    results_path + "associationByLiteratureIndirectOverYears"
)
associationByGeneticIndirectOverYears_file = (
    results_path + "associationByGeneticIndirectOverYears"
)
associationByClinicalIndirectOverYears_file = (
    results_path + "associationByClinicalIndirectOverYears"
)

## plots path
plots_path = results_path + "plots/"


# Prerequired functions
def getDatasourceToName():
    """
    Returns list of data sources weights for overall score.
    """

    names = [[datasource["id"], datasource["label"]] for datasource in dataSources]

    return names


def getDatasourceToColor(palette):
    """
    Returns colourmap dictionary of data sources.
    """

    sources = [name for _, name in getDatasourceToName()]
    colours = {}
    for source, colour in zip(
        sources, sns.color_palette(palette=palette, n_colors=len(sources)).as_hex()
    ):
        colours[source] = colour

    return colours


def getDatatypeToColor():
    """
    Returns colourmap dictionary of data types.
    """

    colours = {
        "all": "#808080",
        "Somatic mutations": "#1f78b4",
        "Literature": "#fdbf6f",
        "Affected pathway": "#cab2d6",
        "Known drug": "#b2df8a",
        "Genetic association": "#fdbf6f",
        "Animal model": "#1f78b4",
        "RNA expression": "#6a3d9a",
    }

    return colours


def plotDiseaseTargetNovelty(
    targetId,
    diseaseId,
    includeNonDated=True,
    major="novelty",
    minor="score",
    showDatasource=True,
    showOverall=True,
    showGenetic=False,
    showClinical=False,
    showLiterature=False,
    img=None,
    vlines=[],
    hlines=[],
    overall_data=associationByOverallIndirectOverYears_file,
    datasource_data=associationByDatasourceIndirectOverYears_file,
    genetic_data=associationByGeneticIndirectOverYears_file,
    clinical_data=associationByClinicalIndirectOverYears_file,
    literature_data=associationByLiteratureIndirectOverYears_file,
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
        overall_data (str):     path to overall score and novelty data
        datasource_data (str):  path to datasource score and novelty data

    Returns:
        img with timeseries plot

    """

    # Establish spark connection
    spark = SparkSession.builder.getOrCreate()

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
        )
    ).toPandas()

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
        .join(datasourceNames, "datasourceId", "left").orderBy("datasourceName")
    ).toPandas()

    # genetic data
    genetic_data = (
        spark.read.parquet(genetic_data)
        # filter disease and target
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
    ).toPandas()

    # literature data
    literature_data = (
        spark.read.parquet(literature_data)
        # filter disease and target
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
    ).toPandas()

    # clinical data
    clinical_data = (
        spark.read.parquet(clinical_data)
        # filter disease and target
        .filter((F.col("targetId") == targetId) & (F.col("diseaseId") == diseaseId))
    ).toPandas()

    # include non-dated evidence labeled as lastYear + 1
    if includeNonDated:
        overall_data.year = overall_data.year.fillna(lastYear + 1)
        datasource_data.year = datasource_data.year.fillna(lastYear + 1)
        clinical_data.year = clinical_data.year.fillna(lastYear + 1)
        literature_data.year = literature_data.year.fillna(lastYear + 1)
        genetic_data.year = genetic_data.year.fillna(lastYear + 1)

    # initialize figure
    plt.figure()
    fig, ax = plt.subplots(figsize=figsize)
    matplotlib.rc("font", **font)
    plt.rcParams["savefig.facecolor"] = "white"

    if showDatasource:
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y=major,
            hue="datasourceName",
            palette=dict(getDatasourceToColor(palette="Paired")),
            lw=thickLineWidth,
            marker="o",
            markersize=0,
            linestyle="-",
            ax=ax,
            legend=True,
        )
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y=minor,
            hue="datasourceName",
            palette=dict(getDatasourceToColor(palette="Paired")),
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )
    if showOverall:
        sns.lineplot(
            data=overall_data.fillna(0),
            x="year",
            y=major,
            color=overallColor,
            lw=thickLineWidth,
            marker="o",
            markersize=0,
            linestyle="-",
            ax=ax,
            label="Overall",
        )
        sns.lineplot(
            data=overall_data.fillna(0),
            x="year",
            y=minor,
            color=overallColor,
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )
    if showClinical:
        sns.lineplot(
            data=clinical_data.fillna(0),
            x="year",
            y=major,
            color=clinicalColor,
            lw=thickLineWidth,
            marker="o",
            markersize=0,
            linestyle="-",
            ax=ax,
            label="Clinical",
        )
        sns.lineplot(
            data=clinical_data.fillna(0),
            x="year",
            y=minor,
            color=clinicalColor,
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )
    if showGenetic:
        sns.lineplot(
            data=genetic_data.fillna(0),
            x="year",
            y=major,
            color=geneticColor,
            lw=thickLineWidth,
            marker="o",
            markersize=0,
            linestyle="-",
            ax=ax,
            label="Genetic",
        )
        sns.lineplot(
            data=genetic_data.fillna(0),
            x="year",
            y=minor,
            color=geneticColor,
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )
    if showLiterature:
        sns.lineplot(
            data=literature_data.fillna(0),
            x="year",
            y=major,
            color=literatureColor,
            lw=thickLineWidth,
            marker="o",
            markersize=0,
            linestyle="-",
            ax=ax,
            label="Literature",
        )
        sns.lineplot(
            data=literature_data.fillna(0),
            x="year",
            y=minor,
            color=literatureColor,
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )

    # legend
    leg = ax.legend(loc="upper left", frameon=False, fontsize=font["size"])
    for line in leg.get_lines():
        line.set_linewidth(thickLineWidth)

    # title
    ax.set_title(
        "Data supporting {} and {} association".format(
            overall_data.diseaseName.values[0].capitalize(),
            overall_data.targetSymbol.values[0],
        ),
        fontsize=font["size"],
    )

    # customize plot
    ax.set_xlabel("")
    if major == "score":
        ylabel = "Association score"
    elif major == "novelty":
        ylabel = "Association novelty"
    ax.set_ylabel(ylabel)
    if includeNonDated:
        ax.set_xlim(firstYear, lastYear + 1)
        ax.text(lastYear - 6, 0.9, "Non-dated:")
    else:
        ax.set_xlim(firstYear, lastYear)
    ax.set_ylim(ymin, ymax)

    # extra features
    for vline in vlines:
        ax.axvline(x=vline, lw=0.5, linestyle="--", color="k")
    for hline in hlines:
        ax.axhline(y=hline, lw=0.5, linestyle="--", color="k")
    if includeNonDated:
        ax.axvline(x=lastYear, lw=0.5, linestyle="-", color="k")

    if img is not None:
        fig.savefig(img, bbox_inches="tight", dpi=300)
        print(img)
        plt.close()

    # return datasource_data.fillna(0)
