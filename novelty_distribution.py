#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "10 Jul 2025"

"""
novelty_distribution.py: Distribution of novel target-disease associations in the Open Targets Platform.
"""

# Setup Google Cloud machine and sumbit job:
# gcloud dataproc clusters create cf-timeseries --image-version 2.2 --region europe-west1 --master-machine-type n1-standard-2 --secondary-worker-type spot --worker-machine-type n1-standard-4 --worker-boot-disk-size 500 --autoscaling-policy=otg-etl --optional-components=JUPYTER --enable-component-gateway --project open-targets-eu-dev
# gcloud dataproc jobs submit pyspark novelty_distribution.py --cluster=cf-timeseries --project=open-targets-eu-dev --region="europe-west1"


from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window, DataFrame
from functools import reduce

# Paths
ot_platform_version = "25.03"

data_path = "gs://open-targets-data-releases/{}/output/".format(ot_platform_version)
disease_file = data_path + "disease"
target_file = data_path + "target"
target_prioritisation_file = data_path + "target_prioritisation"

results_path = "gs://ot-team/cfalaguera/{}/".format(ot_platform_version)
association_by_datasource_dated_file = results_path + "association_by_datasource_dated"

# Functions
prioritized_therapeutic_area = [
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


def get_therapeutic_area_for_disease():
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
            for ranking, therapeuticArea in enumerate(prioritized_therapeutic_area)
        ],
        schema=["therapeuticArea", "ranking"],
    )

    # select top TA for each disease
    partition = Window.partitionBy("diseaseId").orderBy(F.col("ranking"))
    diseases = spark.read.parquet(disease_file).persist()
    diseases = (
        # get diseases mapped to multiple TAs
        diseases.select(
            F.col("id").alias("diseaseId"),
            F.explode("therapeuticAreas").alias("therapeuticArea"),
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


spark = SparkSession.builder.getOrCreate()

# Novels in 2025
if 1:
    associations = (
        spark.read.parquet(association_by_datasource_dated_file)
        .select("targetId", "diseaseId")
        .distinct()
        .join(
            get_therapeutic_area_for_disease(),
            "diseaseId",
            "inner",
        )
        .filter(
            ~F.col("therapeuticArea").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .filter(
            ~F.col("diseaseId").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
    )
    print(
        "associations:",
        associations.count(),
    )
    print(
        "targets:",
        associations.select("targetId").distinct().count(),
    )

    novels = (
        spark.read.parquet(association_by_datasource_dated_file)
        .filter((F.col("novelty") >= 0.1) & (F.col("year") == 2025))
        .select("targetId", "diseaseId")
        .distinct()
        .join(
            get_therapeutic_area_for_disease(),
            "diseaseId",
            "inner",
        )
        .filter(
            ~F.col("therapeuticArea").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .filter(
            ~F.col("diseaseId").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
    )
    print(
        "novel associations in 2025:",
        novels.count(),
    )
    print(
        "novel targets in 2025:",
        novels.select("targetId").distinct().count(),
    )
    print(
        "novel targets in 2025 by biotype:",
        novels.select("targetId")
        .join(
            spark.read.parquet(target_file).select(
                F.col("id").alias("targetId"), "biotype"
            ),
            on="targetId",
            how="left",
        )
        .groupby("biotype")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )
    print(
        "novel targets in 2025 not clinically explored:",
        novels.select("targetId")
        .distinct()
        .join(
            spark.read.parquet(
                association_by_datasource_dated_file + "/sourceId=chembl"
            )
            .select("targetId")
            .distinct(),
            how="anti",
            on="targetId",
        )
        .count(),
    )
    print(
        "novel targets in 2025 with binding ligand:",
        novels.select("targetId")
        .join(
            spark.read.parquet(target_prioritisation_file).select(
                "targetId", "hasLigand"
            ),
            on="targetId",
            how="left",
        )
        .groupby("hasLigand")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )
    print(
        "novel targets in 2025 with safety event:",
        novels.select("targetId")
        .join(
            spark.read.parquet(target_prioritisation_file).select(
                "targetId", "hasSafetyEvent"
            ),
            on="targetId",
            how="left",
        )
        .groupby("hasSafetyEvent")
        .agg(F.size(F.collect_set("targetId")).alias("target"))
        .toPandas(),
    )

    print("novel associations in 2025 by TA:")
    print(
        novels.groupby("therapeuticAreaName")
        .agg(
            F.size(F.collect_set(F.concat("targetId", "diseaseId"))).alias(
                "association"
            )
        )
        .orderBy("association")
        .toPandas(),
    )

    print("novel associations in 2025 by data source:")
    print(
        spark.read.parquet(association_by_datasource_dated_file)
        .filter((F.col("novelty") >= 0.1) & (F.col("year") == 2025))
        .join(
            get_therapeutic_area_for_disease(),
            "diseaseId",
            "inner",
        )
        .filter(
            ~F.col("therapeuticArea").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .filter(
            ~F.col("diseaseId").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .groupby("datasourceId")
        .agg(
            F.size(F.collect_set(F.concat("targetId", "diseaseId"))).alias(
                "association"
            )
        )
        .orderBy("association")
        .toPandas(),
    )

    print("novel associations in 2025 by data source and therapeutic area:")
    print(
        spark.read.parquet(association_by_datasource_dated_file)
        .filter((F.col("novelty") >= 0.1) & (F.col("year") == 2025))
        .join(
            get_therapeutic_area_for_disease(),
            "diseaseId",
            "inner",
        )
        .filter(
            ~F.col("therapeuticArea").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .filter(
            ~F.col("diseaseId").isin(
                [
                    "GO_0008150",  # biological process
                    "EFO_0001444",  # measurement
                    "EFO_0002571",  # medical procedure
                    "EFO_0000651",  # phenotype
                    "EFO_0005932",  # animal disease
                ]
            )
        )
        .groupby("datasourceId", "therapeuticAreaName")
        .agg(
            F.size(F.collect_set(F.concat("targetId", "diseaseId"))).alias(
                "association"
            )
        )
        .orderBy("association")
        .toPandas(),
    )

# Novels distribution
for type in (
    "",
    "_indirect",
):
    for cutoff in [0.1]:
        print(spark.read.parquet(association_by_datasource_dated_file + type).show())
        data = (
            spark.read.parquet(association_by_datasource_dated_file + type)
            .filter(
                (F.col("year") >= 2000)
                & (F.col("year") <= 2025)
                & (F.col("novelty") >= cutoff)
            )
            .join(
                get_therapeutic_area_for_disease(),
                "diseaseId",
                "inner",
            )
            .filter(
                ~F.col("therapeuticArea").isin(
                    [
                        "GO_0008150",  # biological process
                        "EFO_0001444",  # measurement
                        "EFO_0002571",  # medical procedure
                        "EFO_0000651",  # phenotype
                        "EFO_0005932",  # animal disease
                    ]
                )
            )
            .filter(
                ~F.col("diseaseId").isin(
                    [
                        "GO_0008150",  # biological process
                        "EFO_0001444",  # measurement
                        "EFO_0002571",  # medical procedure
                        "EFO_0000651",  # phenotype
                        "EFO_0005932",  # animal disease
                    ]
                )
            )
            .join(
                spark.read.parquet(target_file)
                # .filter(F.col("biotype") == "protein_coding")
                .select(F.col("id").alias("targetId")),
                "targetId",
                "inner",
            )
            .withColumn(
                "therapeuticAreaName",
                F.when(
                    (F.col("therapeuticAreaName") == "cancer or benign tumor"),
                    F.lit("Oncological"),
                )
                # .when(
                #     (
                #         F.col("therapeuticAreaName")
                #         == "genetic, familial or congenital disease"
                #     ),
                #     F.lit("Congenital"),
                # )
                .otherwise(F.lit("Non-oncological")),
            )
            # .withColumn(
            #     "datasourceId",
            #     F.when(
            #         F.col("datasourceId").isin(
            #             [
            #                 "sysbio",
            #                 "crispr",
            #                 "crispr_screen",
            #                 "slapenrich",
            #                 "gene_burden",
            #                 "cancer_biomarkers",
            #                 "reactome",
            #             ]
            #         ),
            #         F.lit("other"),
            #     ).otherwise(F.col("datasourceId")),
            # )
            .withColumn(
                "max_novelty",
                F.max("novelty").over(
                    Window.partitionBy("targetId", "diseaseId", "datasourceId")
                ),
            )
            .filter(F.col("novelty") == F.col("max_novelty"))
            .join(
                spark.createDataFrame(
                    get_datatype_for_datasource(),
                    schema=["datasourceId", "datatypeId"],
                ),
                "datasourceId",
                "left",
            )
        )

        # write by data source
        (  # fill empty years
            spark.createDataFrame(
                data=[[r] for r in range(2000, 2026, 1)],
                # data=[[r] for r in range(1900, 2026, 1)],
                schema=["year"],
            )
            .crossJoin(data.select("datasourceId").distinct())
            .crossJoin(data.select("therapeuticAreaName").distinct())
            .join(
                data.groupby("year", "datasourceId", "therapeuticAreaName").agg(
                    F.size(F.collect_set(F.concat("diseaseId", "targetId"))).alias(
                        "association"
                    )
                ),
                on=["datasourceId", "year", "therapeuticAreaName"],
                how="left",
            )
            .join(
                (
                    data.groupby("datasourceId", "targetId")
                    .agg(F.min("year").alias("year"))
                    .groupby("year", "datasourceId")
                    .agg(F.size(F.collect_set("targetId")).alias("target"))
                ),
                on=["datasourceId", "year"],
                how="left",
            )
            .fillna(0, subset=["association", "target"])
            .write.partitionBy("datasourceId")
            .parquet(
                results_path
                + "novelty_distribution{}_by_datasourceId/cutoff={}".format(
                    type, cutoff
                )
            )
        )

        # write by data type
        (  # fill empty years
            spark.createDataFrame(
                data=[[r] for r in range(2000, 2026, 1)],
                # data=[[r] for r in range(1900, 2026, 1)],
                schema=["year"],
            )
            .crossJoin(data.select("datatypeId").distinct())
            .crossJoin(data.select("therapeuticAreaName").distinct())
            .join(
                data.groupby("year", "datatypeId", "therapeuticAreaName").agg(
                    F.size(F.collect_set(F.concat("diseaseId", "targetId"))).alias(
                        "association"
                    )
                ),
                on=["datatypeId", "year", "therapeuticAreaName"],
                how="left",
            )
            .join(
                (
                    data.groupby("datatypeId", "targetId")
                    .agg(F.min("year").alias("year"))
                    .groupby("year", "datatypeId")
                    .agg(F.size(F.collect_set("targetId")).alias("target"))
                ),
                on=["datatypeId", "year"],
                how="left",
            )
            .fillna(0, subset=["association", "target"])
            .write.partitionBy("datatypeId")
            .parquet(
                results_path
                + "novelty_distribution{}_by_datatypeId/cutoff={}".format(type, cutoff)
            )
        )

        print("cutoff={} and type={}:".format(cutoff, type))
        print()

        print(
            cutoff,
            data.filter(F.col("datasourceId") != "chembl")
            .select("targetId")
            .distinct()
            .count(),
        )
