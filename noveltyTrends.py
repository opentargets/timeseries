#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "15 Oct 2024"

# generate gcloud machine
"""
noveltyTrends.py: Analyse temporal trends in novel drug targets discovery since 2000.
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark noveltyValidation.py --cluster=cf-timeseries2 --project=open-targets-eu-dev --region="europe-west1"
"""

from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window
import itertools
import pyspark
from functools import reduce


# Paths
ot_platform_version = "23.06"
evidence_file = "gs://open-targets-data-releases/{}/output/etl/parquet/evidence".format(
    ot_platform_version
)
evidenceIndirect_file = "gs://ot-team/cfalaguera/novelty/{}/evidenceIndirect".format(
    ot_platform_version
)
molecule_file = "gs://open-targets-data-releases/{}/output/etl/parquet/molecule".format(
    ot_platform_version
)
associationByDatasourceIndirectOverYears_file = "gs://ot-team/cfalaguera/novelty/{}/associationByDatasourceIndirectOverYears".format(
    ot_platform_version
)
associationByDatasourceDirectOverYears_file = (
    "gs://ot-team/cfalaguera/novelty/{}/associationByDatasourceDirectOverYears".format(
        ot_platform_version
    )
)
associationByOverallIndirectOverYears_file = (
    "gs://ot-team/cfalaguera/novelty/{}/associationByOverallIndirectOverYears".format(
        ot_platform_version
    )
)
associationByOverallDirectOverYears_file = (
    "gs://ot-team/cfalaguera/novelty/{}/associationByOverallDirectOverYears".format(
        ot_platform_version
    )
)
evidenceIndirectDated_file = (
    "gs://ot-team/cfalaguera/novelty/{}/evidenceIndirectDated".format(
        ot_platform_version
    )
)
evidenceDirectDated_file = (
    "gs://ot-team/cfalaguera/novelty/{}/evidenceDirectDated".format(ot_platform_version)
)

# Required dataframes
spark = SparkSession.builder.getOrCreate()

prioritizedTherapeuticArea = [
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
        "id": "crispr_manuscript",
        "sectionId": "crisprManuscript",
        "label": "Manuscript",
        "aggregation": "Affected pathway",
        "aggregationId": "affected_pathway",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "",
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


def getDatatypeToName(output="list"):
    id2label_dict = {
        datasource["aggregationId"]: datasource["aggregation"]
        for datasource in dataSources
    }
    id2label_dict["all"] = "All"
    id2label_dict["tcrd"] = "TCRD"
    if output == "dict":
        return id2label_dict
    elif output == "list":
        id2label_list = [[idx, id2label_dict[idx]] for idx in id2label_dict]
        return id2label_list


# novel drug targets peaks
if 1:
    novelDrugTargets = (
        spark.read.option("header", "true")
        .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_34/drugs.csv")
        .join(
            spark.read.parquet(
                "gs://open-targets-data-releases/23.06/output/etl/parquet/evidence/sourceId=chembl"
            )
            .filter(F.col("clinicalPhase") == 4)
            .select("drugId", "targetId", "diseaseId")
            .distinct(),
            "drugId",
            "left",
        )
        .filter(
            F.col("targetId").isNotNull() & (F.col("firstApprovalYear").isNotNull())
        )
        .withColumn(
            "minApprovalYear",
            F.min("firstApprovalYear").over(Window.partitionBy("targetId")),
        )
        .filter(F.col("minApprovalYear") == F.col("firstApprovalYear"))
        .select(
            "drugId",
            "drug",
            "targetId",
            "diseaseId",
            F.col("firstApprovalYear").alias("approvalYear"),
        )
        .distinct()
        # .join(
        #     getTherapeuticAreaForDisease()
        #     .select("diseaseId", "therapeuticAreaName")
        #     .withColumn(
        #         "therapeuticAreaName",
        #         F.when(
        #             F.col("therapeuticAreaName") == "cancer or benign tumor",
        #             "oncology",
        #         ).otherwise("non-oncology"),
        #     ),
        #     "diseaseId",
        #     "left",
        # )
    )

    results = (
        novelDrugTargets
        # source peaks
        .join(
            spark.read.parquet(associationByDatasourceIndirectOverYears_file).select(
                "targetId", "diseaseId", "year", "novelty", "datasourceId"
            ),
            ["diseaseId", "targetId"],
            "left",
        )
        .filter((F.col("novelty") > 0.0) & (F.col("year").isNotNull()))
        .withColumn(
            "maxNovelty",
            F.max("novelty").over(Window.partitionBy("targetId", "datasourceId")),
        )
        .filter(F.col("novelty") == F.col("maxNovelty"))
        .drop("maxNovelty")
    )
    print(
        results.join(
            spark.createDataFrame(
                getDatatypeForDatasource(),
                schema=["datasourceId", "datatypeId"],
            ),
            "datasourceId",
            "left",
        )
        .filter((F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2022))
        .groupBy("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .show()
    )
    print(
        results.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2022)
        )
        .select("targetId")
        .distinct()
        .count()
    )
    results = (
        results.unionByName(
            # clinical phase peaks
            novelDrugTargets.join(
                spark.read.parquet(evidenceIndirectDated_file)
                .filter(F.col("datasourceId") == "chembl")
                .withColumn(
                    "clinicalPhase",
                    F.when(
                        F.col("clinicalPhase").isin(["1.0", "2.0"]),
                        F.lit("1.0/2.0"),
                    ).otherwise(F.col("clinicalPhase")),
                )
                .groupby("targetId", "clinicalPhase", "drugId")
                .agg(F.min("year").alias("year")),
                ["targetId", "drugId"],
                "inner",
            )
            .withColumn("novelty", F.lit(1))
            .withColumnRenamed("clinicalPhase", "datasourceId")
        )
        # final touches
        .withColumn("window", F.col("year") - F.col("approvalYear"))
        .withColumn(
            "discard",
            F.when(
                (F.col("approvalYear") < F.col("year"))
                & (F.col("datasourceId").isin(["1.0/2.0", "3.0"])),
                F.lit(True),
            ).otherwise(F.lit(False)),
        )
        .filter(F.col("discard") == False)
        .join(
            spark.createDataFrame(
                getDatatypeForDatasource(),
                schema=["datasourceId", "datatypeId"],
            ),
            "datasourceId",
            "left",
        )
        .withColumn(
            "datatypeId",
            F.when(F.col("datatypeId").isNotNull(), F.col("datatypeId")).otherwise(
                F.col("datasourceId")
            ),
        )
    )

    print(
        novelDrugTargets.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2022)
        )
        .select("targetId")
        .distinct()
        .count()
    )

    # bell plots
    if 1:
        for minYear, maxYear in [
            (2020, 2022),
            # (2010, 2015),
            # (2008, 2010),
            (2010, 2012),
            # (2000, 2005),
            # (2006, 2010),
            # (2010, 2011),
            # (2012, 2013),
            # (2014, 2015),
            # (2016, 2017),
            # (2018, 2019),
            # (2020, 2021),
            # (2022, 2023),
        ]:

            print(
                minYear,
                maxYear,
                novelDrugTargets.filter(
                    (F.col("approvalYear") >= minYear)
                    & (F.col("approvalYear") <= maxYear)
                )
                .select("targetId")
                .distinct()
                .count(),
            )

            (
                results.replace("somatic_mutation", "somatic_genetic")
                .replace("genetic_association", "somatic_genetic")
                .filter(
                    (F.col("approvalYear") >= minYear)
                    & (F.col("approvalYear") <= maxYear)
                )
                .withColumn(
                    "<-20",
                    F.when(
                        F.col("year") < (F.col("approvalYear") - 20), F.lit(1)
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-20,-16)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 20))
                        & (F.col("year") < (F.col("approvalYear") - 16)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-16,-12)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 16))
                        & (F.col("year") < (F.col("approvalYear") - 12)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-12,-8)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 12))
                        & (F.col("year") < (F.col("approvalYear") - 8)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-8,-4)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 8))
                        & (F.col("year") < (F.col("approvalYear") - 4)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-4,0)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 4))
                        & (F.col("year") < (F.col("approvalYear") - 0)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "0",
                    F.when(
                        F.col("year") == F.col("approvalYear"),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(0,4]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 4))
                        & (F.col("year") > (F.col("approvalYear") + 0)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(4,8]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 8))
                        & (F.col("year") > (F.col("approvalYear") + 4)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(8,12]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 12))
                        & (F.col("year") > (F.col("approvalYear") + 8)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(12,16]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 16))
                        & (F.col("year") > (F.col("approvalYear") + 12)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(16,20]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 20))
                        & (F.col("year") > (F.col("approvalYear") + 16)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    ">20",
                    F.when(
                        F.col("year") > F.col("approvalYear") + 20, F.lit(1)
                    ).otherwise(F.lit(0)),
                )
                # .withColumn("anytime", F.lit(1))
                .groupby("targetId", "datatypeId")
                .agg(
                    F.max("<-20").alias("<-20"),
                    F.max("[-20,-16)").alias("[-20,-16)"),
                    F.max("[-16,-12)").alias("[-16,-12)"),
                    F.max("[-12,-8)").alias("[-12,-8)"),
                    F.max("[-8,-4)").alias("[-8,-4)"),
                    F.max("[-4,0)").alias("[-4,0)"),
                    F.max("0").alias("0"),
                    F.max("(0,4]").alias("(0,4]"),
                    F.max("(4,8]").alias("(4,8]"),
                    F.max("(8,12]").alias("(8,12]"),
                    F.max("(12,16]").alias("(12,16]"),
                    F.max("(16,20]").alias("(16,20]"),
                    F.max(">20").alias(">20"),
                    # F.max("anytime").alias("anytime"),
                )
                .groupby("datatypeId")
                .sum()
                .toPandas()
                .set_index("datatypeId")
                / (
                    novelDrugTargets.filter(
                        (F.col("approvalYear") >= minYear)
                        & (F.col("approvalYear") <= maxYear)
                    )
                    .select("targetId")
                    .distinct()
                    .count()
                )
                * 100
            ).reset_index().melt(
                id_vars=["datatypeId"], var_name="years", value_name="%withPeak"
            ).to_csv(
                "gs://ot-team/cfalaguera/noveltyBenchmark/bells/years={}-{}".format(
                    minYear,
                    maxYear,
                )
            )

    # cascades plots
    if 1:
        for minYear, maxYear in [
            (2000, 2005),
            (2006, 2010),
            (2011, 2012),
            (2013, 2014),
            (2015, 2016),
            (2017, 2018),
            (2019, 2020),
            (2021, 2022),
        ]:

            results.filter(
                (F.col("approvalYear") >= minYear) & (F.col("approvalYear") <= maxYear)
            ).write.parquet(
                "gs://ot-team/cfalaguera/noveltyBenchmark/cascadesBoxplot/years={}-{}".format(
                    minYear,
                    maxYear,
                )
            )
            print(
                minYear,
                maxYear,
                results.filter(
                    (F.col("approvalYear") >= minYear)
                    & (F.col("approvalYear") <= maxYear)
                )
                .select("targetId")
                .distinct()
                .count(),
            )
