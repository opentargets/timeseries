#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "24 May 2023"

"""
plotTimeries.py: Plot timeseries for associations in Open Targets Platform.
"""

import matplotlib
from matplotlib import pyplot as plt
from pyspark.sql import functions as F
from pyspark.sql import SparkSession
import seaborn as sns


# Settings
font = {"weight": "normal", "size": 15}
overallColor = "lightgray"
firstYear = 1960
lastYear = 2023
ymin = 0
ymax = 1
figsize = (15, 4)
cutoffLinewidth = 0.1
thickLineWidth = 2.5
slimLineWidth = 0.5

# Paths
ot_platform_version = "/23.02/"

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

## plots path
plots_path = results_path + "plots/"


# Prerequired functions
def getDatasourceToName():
    """
    Returns list of data sources weights for overall score.
    """

    weights = [
        ["cancer_biomarkers", "Cancer Biomarkers"],
        ["cancer_gene_census", "Cancer Gene Census"],
        ["chembl", "ChEMBL"],
        ["clingen", "Clingen"],
        ["crispr", "Project Score"],
        ["encore", "ENCORE"],
        ["europepmc", "Europe PMC"],
        ["eva", "ClinVar"],
        ["eva_somatic", "ClinVar (somatic)"],
        ["expression_atlas", "Expression Atlas"],
        ["gene2phenotype", "Gene2phenotype"],
        ["gene_burden", "Gene Burden"],
        ["genomics_england", "GEL PanelApp"],
        ["impc", "IMPC"],
        ["intogen", "IntOGen"],
        ["orphanet", "Orphanet"],
        ["ot_crispr", "OT CRISPR"],
        ["ot_crispr_validation", "OT Validation"],
        ["ot_genetics_portal", "OT Genetics"],
        ["progeny", "PROGENy"],
        ["reactome", "Reactome"],
        ["slapenrich", "SLAPenrich"],
        ["sysbio", "Gene signatures"],
        ["uniprot_literature", "UniProt literature"],
        ["uniprot_variants", "UniProt curated variants"],
    ]

    return weights


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


def plotDiseaseTargetNovelty(
    targetId,
    diseaseId,
    score=False,
    novelty=True,
    img=None,
    label=0,
    vlines=[],
    hlines=[],
    overall_data=associationByOverallIndirectOverYears_file,
    datasource_data=associationByDatasourceIndirectOverYears_file,
):
    """
    Plot timeseries for a disease-target associations.

    Args:
        targetId (str):         target
        diseaseId (str):        disease
        score (bool):           plot association score evolution
        novelty (bool):         plot association novelty evolution
        img (str):              path to save image
        label (str):            label to write on the plot
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

    # initialize figure
    plt.figure()
    fig, ax = plt.subplots(figsize=figsize)
    matplotlib.rc("font", **font)
    plt.rcParams["savefig.facecolor"] = "white"

    # plot datasource score
    if score:
        sns.lineplot(
            data=datasource_data.fillna(0),
            x="year",
            y="score",
            hue="datasourceName",
            palette=dict(getDatasourceToColor(palette="Paired")),
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )

    # plot datasource novelty
    sns.lineplot(
        data=datasource_data.fillna(0),
        x="year",
        y="novelty",
        hue="datasourceName",
        palette=dict(getDatasourceToColor(palette="Paired")),
        lw=thickLineWidth,
        marker="o",
        markersize=0,
        linestyle="-",
        ax=ax,
        legend=True,
    )

    # plot overall score
    if score:
        sns.lineplot(
            data=overall_data.fillna(0),
            x="year",
            y="score",
            color=overallColor,
            lw=slimLineWidth,
            marker="o",
            markersize=0,
            linestyle="--",
            ax=ax,
            legend=False,
        )

    # plot overall novelty
    sns.lineplot(
        data=overall_data.fillna(0),
        x="year",
        y="novelty",
        color=overallColor,
        lw=thickLineWidth,
        marker="o",
        markersize=0,
        linestyle="-",
        ax=ax,
        legend=True,
        label="Overall",
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
    ax.set_ylabel("Novel data")
    ax.set_xlim(firstYear, lastYear)
    ax.set_ylim(ymin, ymax)

    # extra features
    for vline in vlines:
        ax.axvline(x=vline, lw=0.5, linestyle="--", color="k")
    for hline in hlines:
        ax.axhline(y=hline, lw=0.5, linestyle="--", color="k")
    if label:
        ax.text(1951, 0.1, label)

    if img is not None:
        fig.savefig(img, bbox_inches="tight", dpi=300)
        print(img)
        plt.close()

    # return datasource_data.fillna(0)
    return 0


def getDiseaseTargetSignature(diseaseId, targetId):
    """
    Get novelty signatures for a disease-target associations.

    Args:
        targetId (str):         target
        diseaseId (str):        disease

    Returns:
        print novelty signatures by datasource
    """

    # Establish spark connection
    spark = SparkSession.builder.getOrCreate()

    # datasource names
    datasourceNames = spark.createDataFrame(
        data=[
            [datasourceId, datasourceName]
            for datasourceId, datasourceName in getDatasourceToName()
        ],
        schema=["datasourceId", "datasourceName"],
    )

    # print signatures
    for datasourceName, noveltySignature in (
        spark.read.parquet(associationByDatasourceIndirectOverYearsSignature_file)
        # filter disease and target
        .filter((F.col("diseaseId") == diseaseId) & (F.col("targetId") == targetId))
        # add datasources' fancy names
        .join(datasourceNames, "datasourceId", "left")
        .orderBy("datasourceName")
        .select("datasourceName", "noveltySignature")
        .toPandas()
        .values
    ):
        print("".join(map(str, datasourceName)))
        print("".join(map(str, noveltySignature)))
