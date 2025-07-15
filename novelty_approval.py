#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "15 Jul 2025"

# generate gcloud machine
"""
novelty_approval.py: Analyse temporal trends in novel drug targets discovery since 2000.
"""

# Setup Google Cloud machine and sumbit job:
# gcloud dataproc clusters create cf-timeseries --image-version 2.2 --region europe-west1 --master-machine-type n1-standard-2 --secondary-worker-type spot --worker-machine-type n1-standard-4 --worker-boot-disk-size 500 --autoscaling-policy=otg-etl --optional-components=JUPYTER --enable-component-gateway --project open-targets-eu-dev
# gcloud dataproc jobs submit pyspark novelty_approval.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"

from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window
from pyspark.sql import types as T


# Paths
ot_platform_version = "25.03"
association_by_datasource_dated_indirect_file = (
    "gs://ot-team/cfalaguera/{}/association_by_datasource_dated_indirect".format(
        ot_platform_version
    )
)
association_by_datasource_dated_file = (
    "gs://ot-team/cfalaguera/{}/association_by_datasource_dated".format(
        ot_platform_version
    )
)
evidence_dated_indirect_file = (
    "gs://ot-team/cfalaguera/{}/evidence_dated_indirect".format(ot_platform_version)
)
evidence_dated_file = "gs://ot-team/cfalaguera/{}/evidence_dated".format(
    ot_platform_version
)
chembl_evidence = (
    "gs://open-targets-data-releases/25.03/output/evidence/sourceId=chembl"
)

novel_drug_target_file = "gs://ot-team/cfalaguera/{}/novel_drug_target".format(
    ot_platform_version
)


# Required dataframes
spark = SparkSession.builder.getOrCreate()

dataSources = [
    {
        "id": "gwas_credible_sets",
        "sectionId": "gwasCredibleSets",
        "label": "GWAS associations",
        "aggregation": "Genetic association",
        "aggregationId": "genetic_association",
        "weight": 1.0,  # needs to be a float
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#gwas-credible-sets",
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


def get_datatype_for_datasource(output="list"):
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


# Novelty profiles for drug targets since 2000 (indirect)
if 1:
    novelDrugTargets = (
        spark.read.option("header", "true")
        # .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_34/drugs.csv")
        .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_35/drugs.csv")
        .join(
            spark.read.parquet(chembl_evidence)
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
    )

    if 0:
        novelDrugTargets.select("targetId", "diseaseId").distinct().write.parquet(
            novel_drug_target_file
        )

    results = (
        novelDrugTargets
        # source peaks
        .join(
            spark.read.parquet(association_by_datasource_dated_indirect_file).select(
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
                get_datatype_for_datasource(),
                schema=["datasourceId", "datatypeId"],
            ),
            "datasourceId",
            "left",
        )
        .filter((F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025))
        .groupBy("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .show()
    )
    print(
        results.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025)
        )
        .select("targetId")
        .distinct()
        .count()
    )
    results = (
        results.unionByName(
            # clinical phase peaks
            novelDrugTargets.join(
                spark.read.parquet(evidence_dated_indirect_file)
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
                get_datatype_for_datasource(),
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
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025)
        )
        .select("targetId")
        .distinct()
        .count()
    )

    print(
        "target supported",
        results.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025)
        )
        .withColumn(
            "datatypeId",
            F.when(
                F.col("datatypeId").isin(["genetic_association", "somatic_mutation"]),
                F.lit("genetic"),
            )
            .when(
                F.col("datatypeId").isin(
                    ["affected_pathway", "rna_expression", "animal_model"]
                ),
                F.lit("other"),
            )
            .when(
                F.col("datatypeId") == "literature",
                F.lit("literature"),
            )
            .otherwise(F.col("datatypeId")),
        )
        .groupby("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )

    # bell plots every 4 years
    if 0:
        for minYear, maxYear in [
            # (2020, 2022),
            (2010, 2015),
            # (2008, 2010),
            # (2010, 2012),
            # (2000, 2005),
            # (2006, 2010),
            # (2010, 2011),
            # (2012, 2013),
            # (2014, 2015),
            # (2016, 2017),
            # (2018, 2019),
            # (2020, 2021),
            # (2022, 2023),
            (2000, 2005),
            # (2010, 2015),
            (2020, 2025),
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

            print(results.select("datasourceId").distinct())

            (
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
                )
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
                "gs://ot-team/cfalaguera/{}/bell_indirect/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # bell plots every 5 years
    if 0:
        for minYear, maxYear in [
            # (2020, 2022),
            (2010, 2015),
            # (2008, 2010),
            # (2010, 2012),
            # (2000, 2005),
            # (2006, 2010),
            # (2010, 2011),
            # (2012, 2013),
            # (2014, 2015),
            # (2016, 2017),
            # (2018, 2019),
            # (2020, 2021),
            # (2022, 2023),
            (2000, 2005),
            # (2010, 2015),
            (2020, 2025),
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
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
                )
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
                    "[-20,-15)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 20))
                        & (F.col("year") < (F.col("approvalYear") - 15)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-15,-10)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 15))
                        & (F.col("year") < (F.col("approvalYear") - 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-10,-5)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 10))
                        & (F.col("year") < (F.col("approvalYear") - 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-5,0)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 5))
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
                    "(0,5]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 5))
                        & (F.col("year") > (F.col("approvalYear") + 0)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(5,10]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 10))
                        & (F.col("year") > (F.col("approvalYear") + 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(10,15]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 15))
                        & (F.col("year") > (F.col("approvalYear") + 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(15,20]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 20))
                        & (F.col("year") > (F.col("approvalYear") + 15)),
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
                    F.max("[-20,-15)").alias("[-20,-15)"),
                    F.max("[-15,-10)").alias("[-15,-10)"),
                    F.max("[-10,-5)").alias("[-10,-5)"),
                    F.max("[-5,0)").alias("[-5,0)"),
                    F.max("0").alias("0"),
                    F.max("(0,5]").alias("(0,5]"),
                    F.max("(5,10]").alias("(5,10]"),
                    F.max("(10,15]").alias("(10,15]"),
                    F.max("(15,20]").alias("(15,20]"),
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
                "gs://ot-team/cfalaguera/{}/bell_indirect/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # cascades plots intervals
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2006, 2010),
            (2011, 2012),
            (2013, 2014),
            (2015, 2016),
            (2017, 2018),
            (2019, 2020),
            (2021, 2022),
            (2023, 2024),
            (2025, 2026),
        ]:

            results.filter(
                (F.col("approvalYear") >= minYear) & (F.col("approvalYear") <= maxYear)
            ).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascadesBoxplot/years={}-{}".format(
                    ot_platform_version,
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

    # cascades plots
    if 0:
        for year in range(2000, 2005):
            print(year)

            results.filter((F.col("approvalYear") == year)).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascade/approvalYear={}".format(
                    ot_platform_version, year
                )
            )
            print(
                year,
                results.filter((F.col("approvalYear") == year))
                .select("targetId")
                .distinct()
                .count(),
            )

# Novelty profiles for drug targets since 2000 (direct)
if 1:
    novelDrugTargets = (
        spark.read.option("header", "true")
        # .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_34/drugs.csv")
        .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_35/drugs.csv")
        .join(
            spark.read.parquet(chembl_evidence)
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
        .filter(
            (F.col("firstApprovalYear") >= 2000) & (F.col("firstApprovalYear") <= 2025)
        )
        .select(
            "drugId",
            "drug",
            "targetId",
            "diseaseId",
            F.col("firstApprovalYear").alias("approvalYear"),
        )
        .distinct()
    )

    results = (
        novelDrugTargets
        # source peaks
        .join(
            spark.read.parquet(association_by_datasource_dated_file).select(
                "targetId", "diseaseId", "year", "novelty", "datasourceId"
            ),
            ["diseaseId", "targetId"],
            "left",
        )
        .filter((F.col("novelty") > 0) & (F.col("year").isNotNull()))
        .withColumn(
            "maxNovelty",
            F.max("novelty").over(Window.partitionBy("targetId", "datasourceId")),
        )
        .filter(F.col("novelty") == F.col("maxNovelty"))
        .drop("maxNovelty")
    )

    # print(
    #     results.join(
    #         spark.createDataFrame(
    #             get_datatype_for_datasource(),
    #             schema=["datasourceId", "datatypeId"],
    #         ),
    #         "datasourceId",
    #         "left",
    #     )
    #     .filter((F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025))
    #     .groupBy("datatypeId")
    #     .agg(F.size(F.collect_set("targetId")).alias("targets"))
    #     .show()
    # )

    results = (
        results.unionByName(
            # clinical phase peaks
            novelDrugTargets.join(
                spark.read.parquet(evidence_dated_file)
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
                get_datatype_for_datasource(),
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
        "target",
        novelDrugTargets.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025)
        )
        .select("targetId")
        .distinct()
        .count(),
    )
    print(
        "target supported",
        results.filter(
            (F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025)
        )
        .withColumn(
            "datatypeId",
            F.when(
                F.col("datatypeId").isin(["genetic_association", "somatic_mutation"]),
                F.lit("genetic"),
            )
            .when(
                F.col("datatypeId").isin(
                    ["affected_pathway", "rna_expression", "animal_model"]
                ),
                F.lit("other"),
            )
            .when(
                F.col("datatypeId") == "literature",
                F.lit("literature"),
            )
            .otherwise(F.col("datatypeId")),
        )
        .groupby("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )

    # bell plots every 4 years
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2005, 2010),
            (2010, 2015),
            (2015, 2020),
            (2020, 2025),
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
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
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
                "gs://ot-team/cfalaguera/{}/bell/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # bell plots every 5 years
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2005, 2010),
            (2010, 2015),
            (2015, 2020),
            (2020, 2025),
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
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
                )
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
                    "[-20,-15)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 20))
                        & (F.col("year") < (F.col("approvalYear") - 15)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-15,-10)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 15))
                        & (F.col("year") < (F.col("approvalYear") - 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-10,-5)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 10))
                        & (F.col("year") < (F.col("approvalYear") - 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-5,0)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 5))
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
                    "(0,5]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 5))
                        & (F.col("year") > (F.col("approvalYear") + 0)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(5,10]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 10))
                        & (F.col("year") > (F.col("approvalYear") + 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(10,15]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 15))
                        & (F.col("year") > (F.col("approvalYear") + 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(15,20]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 20))
                        & (F.col("year") > (F.col("approvalYear") + 15)),
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
                    F.max("[-20,-15)").alias("[-20,-15)"),
                    F.max("[-15,-10)").alias("[-15,-10)"),
                    F.max("[-10,-5)").alias("[-10,-5)"),
                    F.max("[-5,0)").alias("[-5,0)"),
                    F.max("0").alias("0"),
                    F.max("(0,5]").alias("(0,5]"),
                    F.max("(5,10]").alias("(5,10]"),
                    F.max("(10,15]").alias("(10,15]"),
                    F.max("(15,20]").alias("(15,20]"),
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
                "gs://ot-team/cfalaguera/{}/bell/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # cascades plots intervals
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2006, 2010),
            (2011, 2012),
            (2013, 2014),
            (2015, 2016),
            (2017, 2018),
            (2019, 2020),
            (2021, 2022),
            (2023, 2024),
            (2025, 2026),
        ]:

            results.filter(
                (F.col("approvalYear") >= minYear) & (F.col("approvalYear") <= maxYear)
            ).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascadesBoxplot/years={}-{}".format(
                    ot_platform_version,
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

    # cascades plots
    if 0:
        for year in range(2000, 2026):
            print(year)

            results.filter((F.col("approvalYear") == year)).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascade_direct/approvalYear={}".format(
                    ot_platform_version, year
                )
            )
            print(
                year,
                results.filter((F.col("approvalYear") == year))
                .select("targetId")
                .distinct()
                .count(),
            )

# Novelty profiles for drug targets since 2000 base on early clinical phase <3 (direct)
if 1:
    novelDrugTargets = (
        spark.read.option("header", "true")
        # .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_34/drugs.csv")
        .csv("gs://ot-team/cfalaguera/ChEMBL/CHEMBL_35/drugs.csv")
        .join(
            spark.read.parquet(chembl_evidence)
            .filter(F.col("clinicalPhase") <= 2.0)
            .select(
                "drugId",
                "targetId",
                "diseaseId",
                F.split(F.col("studyStartDate"), "-")
                .getItem(0)
                .cast(T.IntegerType())
                .alias("clinicalYear"),
            )
            .distinct(),
            "drugId",
            "left",
        )
        .filter(
            F.col("targetId").isNotNull() & (F.col("firstApprovalYear").isNotNull())
        )
        .withColumn(
            "minClinicalYear",
            F.min("clinicalYear").over(Window.partitionBy("targetId")),
        )
        .filter(F.col("minClinicalYear") == F.col("clinicalYear"))
        .filter((F.col("clinicalYear") >= 2000) & (F.col("clinicalYear") <= 2025))
        .select(
            "drugId",
            "drug",
            "targetId",
            "diseaseId",
            F.col("minClinicalYear").alias("clinicalYear"),
        )
        .distinct()
    )

    results = (
        novelDrugTargets
        # source peaks
        .join(
            spark.read.parquet(association_by_datasource_dated_file).select(
                "targetId", "diseaseId", "year", "novelty", "datasourceId"
            ),
            ["diseaseId", "targetId"],
            "left",
        )
        .filter((F.col("novelty") > 0) & (F.col("year").isNotNull()))
        .withColumn(
            "maxNovelty",
            F.max("novelty").over(Window.partitionBy("targetId", "datasourceId")),
        )
        .filter(F.col("novelty") == F.col("maxNovelty"))
        .drop("maxNovelty")
    )

    # print(
    #     results.join(
    #         spark.createDataFrame(
    #             get_datatype_for_datasource(),
    #             schema=["datasourceId", "datatypeId"],
    #         ),
    #         "datasourceId",
    #         "left",
    #     )
    #     .filter((F.col("approvalYear") >= 2000) & (F.col("approvalYear") <= 2025))
    #     .groupBy("datatypeId")
    #     .agg(F.size(F.collect_set("targetId")).alias("targets"))
    #     .show()
    # )

    results = (
        results.unionByName(
            # clinical phase peaks
            novelDrugTargets.join(
                spark.read.parquet(evidence_dated_file)
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
        .withColumn("window", F.col("year") - F.col("clinicalYear"))
        .withColumn(
            "discard",
            F.when(
                (F.col("clinicalYear") < F.col("year"))
                & (F.col("datasourceId").isin(["1.0/2.0", "3.0"])),
                F.lit(True),
            ).otherwise(F.lit(False)),
        )
        .filter(F.col("discard") == False)
        .join(
            spark.createDataFrame(
                get_datatype_for_datasource(),
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
        "target",
        novelDrugTargets.filter(
            (F.col("clinicalYear") >= 2000) & (F.col("clinicalYear") <= 2025)
        )
        .select("targetId")
        .distinct()
        .count(),
    )
    print(
        "target supported",
        results.filter(
            (F.col("clinicalYear") >= 2000) & (F.col("clinicalYear") <= 2025)
        )
        .withColumn(
            "datatypeId",
            F.when(
                F.col("datatypeId").isin(["genetic_association", "somatic_mutation"]),
                F.lit("genetic"),
            )
            .when(
                F.col("datatypeId").isin(
                    ["affected_pathway", "rna_expression", "animal_model"]
                ),
                F.lit("other"),
            )
            .when(
                F.col("datatypeId") == "literature",
                F.lit("literature"),
            )
            .otherwise(F.col("datatypeId")),
        )
        .groupby("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )

    # bell plots every 4 years
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2005, 2010),
            (2010, 2015),
            (2015, 2020),
            (2020, 2025),
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
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
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
                "gs://ot-team/cfalaguera/{}/bell/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # bell plots every 5 years
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2005, 2010),
            (2010, 2015),
            (2015, 2020),
            (2020, 2025),
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
                results.withColumn(
                    "datatypeId",
                    F.when(
                        (F.col("datatypeId") == "somatic_mutation")
                        | (F.col("datatypeId") == "genetic_association"),
                        F.lit("somatic_genetic"),
                    )
                    .when(
                        (F.col("datatypeId") == "affected_pathway")
                        | (F.col("datatypeId") == "animal_model")
                        | (F.col("datatypeId") == "rna_expression"),
                        F.lit("other"),
                    )
                    .when(
                        (F.col("datatypeId") == "literature"),
                        F.col("datatypeId"),
                    )
                    .otherwise(F.col("datasourceId")),
                )
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
                    "[-20,-15)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 20))
                        & (F.col("year") < (F.col("approvalYear") - 15)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-15,-10)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 15))
                        & (F.col("year") < (F.col("approvalYear") - 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-10,-5)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 10))
                        & (F.col("year") < (F.col("approvalYear") - 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "[-5,0)",
                    F.when(
                        (F.col("year") >= (F.col("approvalYear") - 5))
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
                    "(0,5]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 5))
                        & (F.col("year") > (F.col("approvalYear") + 0)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(5,10]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 10))
                        & (F.col("year") > (F.col("approvalYear") + 5)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(10,15]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 15))
                        & (F.col("year") > (F.col("approvalYear") + 10)),
                        F.lit(1),
                    ).otherwise(F.lit(0)),
                )
                .withColumn(
                    "(15,20]",
                    F.when(
                        (F.col("year") <= (F.col("approvalYear") + 20))
                        & (F.col("year") > (F.col("approvalYear") + 15)),
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
                    F.max("[-20,-15)").alias("[-20,-15)"),
                    F.max("[-15,-10)").alias("[-15,-10)"),
                    F.max("[-10,-5)").alias("[-10,-5)"),
                    F.max("[-5,0)").alias("[-5,0)"),
                    F.max("0").alias("0"),
                    F.max("(0,5]").alias("(0,5]"),
                    F.max("(5,10]").alias("(5,10]"),
                    F.max("(10,15]").alias("(10,15]"),
                    F.max("(15,20]").alias("(15,20]"),
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
                "gs://ot-team/cfalaguera/{}/bell/years={}-{}".format(
                    ot_platform_version,
                    minYear,
                    maxYear,
                )
            )

    # cascades plots intervals
    if 0:
        for minYear, maxYear in [
            (2000, 2005),
            (2006, 2010),
            (2011, 2012),
            (2013, 2014),
            (2015, 2016),
            (2017, 2018),
            (2019, 2020),
            (2021, 2022),
            (2023, 2024),
            (2025, 2026),
        ]:

            results.filter(
                (F.col("approvalYear") >= minYear) & (F.col("approvalYear") <= maxYear)
            ).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascadesBoxplot/years={}-{}".format(
                    ot_platform_version,
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

    # cascades plots
    if 1:
        for year in range(2000, 2026):
            print(year)

            results.filter((F.col("clinicalYear") == year)).write.parquet(
                "gs://ot-team/cfalaguera/{}/cascade_direct_by_clinical/clinicalYear={}".format(
                    ot_platform_version, year
                )
            )
