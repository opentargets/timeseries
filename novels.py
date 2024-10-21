#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "15 Oct 2024"

"""
novels.py: List novel target-disease associations in the Open Targets Platform.
"""

# sumbit job to gcloud machine
"""
gcloud dataproc jobs submit pyspark novels.py --cluster=cf-novelty --project=open-targets-eu-dev --region="europe-west1"
"""
from pyspark.sql import functions as F
from pyspark.sql import SparkSession, Window, DataFrame
from functools import reduce

# Required data
year = 2023
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
excludeTherapeuticArea = [
    "GO_0008150",  # biological process
    "EFO_0001444",  # measurement
    "EFO_0002571",  # medical procedure
    # "EFO_0000651",  # phenotype
]
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTchembl",
    },
    {
        "id": "phase=IV",
        "sectionId": "phase=IV",
        "label": "phase=IV",
        "aggregation": "phase=IV",
        "aggregationId": "phase=IV",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#chembl",
        "category": "OTphase=IV",
    },
    {
        "id": "bioactivity",
        "sectionId": "bioactivity",
        "label": "Bioactivity",
        "aggregation": "bioactivity",
        "aggregationId": "bioactivity",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#chembl",
        "category": "bioactivity",
    },
    {
        "id": "phase<IV",
        "sectionId": "phase<IV",
        "label": "phase<IV",
        "aggregation": "phase<IV",
        "aggregationId": "phase<IV",
        "weight": 1.0,
        "isPrivate": False,
        "docsLink": "https://platform-docs.opentargets.org/evidence#chembl",
        "category": "OTphase<IV",
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
        "category": "OTgen",
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
        "category": "OTgen",
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
        "category": "OTpath",
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
        "category": "OTpath",
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
        "category": "OTpath",
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
        "category": "OTpath",
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
        "category": "OTlit",
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
        "category": "OTrna",
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
        "category": "OTgen",
    },
    # {
    #     "id": "ot_crispr",
    #     "sectionId": "oOTgen",
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
prioritisationColumns = [
    {
        "id": "maxClinicalTrialPhase",
        "label": "Target in clinic",
        "aggregation": "aggregations.precedence",
        "sectionId": "knownDrugs",
        "description": "Target is in clinical trials for any indication",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#target-in-clinic",
    },
    {
        "id": "isInMembrane",
        "label": "Membrane protein",
        "aggregation": "aggregations.tractability",
        "sectionId": "subcellularLocation",
        "description": "Target is annotated to be located in the cell membrane",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#membrane-protein",
    },
    {
        "id": "isSecreted",
        "label": "Secreted protein",
        "aggregation": "aggregations.tractability",
        "sectionId": "subcellularLocation",
        "description": "Target is annotated to be secreted",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#secreted-protein",
    },
    {
        "id": "hasLigand",
        "label": "Ligand binder",
        "aggregation": "aggregations.tractability",
        "sectionId": "tractability",
        "description": "Target binds a specific ligand",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#ligand-binder",
    },
    {
        "id": "hasSmallMoleculeBinder",
        "label": "Small molecule binder",
        "aggregation": "aggregations.tractability",
        "sectionId": "tractability",
        "description": "Target binds a small molecule",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#small-molecule-binder",
    },
    {
        "id": "hasPocket",
        "label": "Predicted pockets",
        "aggregation": "aggregations.tractability",
        "sectionId": "tractability",
        "description": "Target has predicted pockets",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#predicted-pockets",
    },
    {
        "id": "mouseOrthologMaxIdentityPercentage",
        "label": "Mouse ortholog identity",
        "aggregation": "aggregations.doability",
        "sectionId": "compGenomics",
        "description": "Mouse ortholog maximum identity percentage",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#mouse-ortholog-identity",
    },
    {
        "id": "hasHighQualityChemicalProbes",
        "label": "Chemical probes",
        "aggregation": "aggregations.doability",
        "sectionId": "chemicalProbes",
        "description": "Availability of high quality chemical probes for the target",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#chemical-probes",
    },
    {
        "id": "mouseKOScore",
        "label": "Mouse models",
        "aggregation": "aggregations.safety",
        "sectionId": "mousePhenotypes",
        "description": "Availability of mouse knockout models for the target",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#mouse-models",
    },
    {
        "id": "geneticConstraint",
        "label": "Genetic constraint",
        "aggregation": "aggregations.safety",
        "sectionId": "geneticConstraint",
        "description": "Relative genetic constraint in natural populations derived from GnomAD",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#genetic-constraint",
    },
    {
        "id": "geneEssentiality",
        "label": "Gene essentiality",
        "aggregation": "aggregations.safety",
        "sectionId": "depMapEssentiality",
        "description": "Gene is defined as core essential by the DepMap portal",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#gene-essentiality",
    },
    {
        "id": "hasSafetyEvent",
        "label": "Known adverse events",
        "aggregation": "aggregations.safety",
        "sectionId": "safety",
        "description": "Target associated with a curated adverse event",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#known-adverse-events",
    },
    {
        "id": "isCancerDriverGene",
        "label": "Cancer driver gene",
        "aggregation": "aggregations.safety",
        "sectionId": "cancerHallmarks",
        "description": "Target is classified as an Oncogene and/or Tumor Suppressor Gene",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#cancer-driver-gene",
    },
    {
        "id": "paralogMaxIdentityPercentage",
        "label": "Paralogues",
        "aggregation": "aggregations.safety",
        "sectionId": "compGenomics",
        "description": "Paralog maximum identity percentage",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#paralogues",
    },
    {
        "id": "tissueSpecificity",
        "label": "Tissue specificity",
        "aggregation": "aggregations.safety",
        "sectionId": "expressions",
        "description": "HPA category types of elevated expression across tissues for the target",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#tissue-specificity",
    },
    {
        "id": "tissueDistribution",
        "label": "Tissue distribution",
        "aggregation": "aggregations.safety",
        "sectionId": "expressions",
        "description": "HPA category types of detectable expression across tissues for the target",
        "docsLink": "https://platform-docs.opentargets.org/target-prioritisation#tissue-distribution",
    },
]

examples = [
    ("ENSG00000237380", "EFO_0003060"),
    # ("ENSG00000149506", "EFO_0008560"),  # bg-
    # ("ENSG00000211891", "MONDO_0005271"),  # bg-
    # ("ENSG00000104687", "EFO_0000519"),  # bg-
    # ("ENSG00000139687", "EFO_0000702"),  # bg-
    # ("ENSG00000141736", "EFO_0003060"),  # bg-
    ("ENSG00000214842", "EFO_0000508"),
    ("ENSG00000205038", "EFO_0008560"),
    ("ENSG00000204544", "EFO_0000519"),
    ("ENSG00000185842", "EFO_0010642"),
    ("ENSG00000184956", "EFO_0000702"),
    ("ENSG00000183578", "EFO_0000181"),
    ("ENSG00000180483", "EFO_0000508"),
    ("ENSG00000177238", "MONDO_0007254"),
    ("ENSG00000165972", "EFO_0004248"),
    ("ENSG00000164675", "EFO_0004248"),
    ("ENSG00000162779", "EFO_0000279"),
    ("ENSG00000154007", "EFO_0000279"),
    ("ENSG00000146453", "EFO_0000279"),
    ("ENSG00000145569", "MONDO_0005271"),
    ("ENSG00000121211", "MONDO_0007254"),
    ("ENSG00000089723", "EFO_0000707"),
    ("ENSG00000064225", "EFO_0003060"),
]

spark = SparkSession.builder.getOrCreate()


def getDatasourceToName(output="list"):
    """
    Returns list of data sources names/labels.
    """

    id2label_dict = {
        datasource["id"]: datasource["label"] for datasource in dataSources
    }
    if output == "dict":
        return id2label_dict
    elif output == "list":
        id2label_list = [[idx, id2label_dict[idx]] for idx in id2label_dict]
        return id2label_list


def getPrioritisationToAggregation(output="list"):

    id2label_dict = {col["id"]: col["aggregation"] for col in prioritisationColumns}
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


# Paths
ot_platform_version = "23.06"
data_path = "gs://open-targets-data-releases/{}/output/etl/parquet/".format(
    ot_platform_version
)
diseases_file = data_path + "diseases"
targets_file = data_path + "targets"
drugs_file = data_path + "searchDrug"
targetsPrioritisation_file = (
    "gs://open-targets-data-releases/24.06/output/etl/parquet/targetPrioritisation"
)
results_path = "gs://ot-team/cfalaguera/novelty/{}/".format(ot_platform_version)
associationByDatasourceIndirectOverYears_file = (
    results_path + "associationByDatasourceIndirectOverYears"
)
associationByDatasourceDirectOverYears_file = (
    results_path + "associationByDatasourceDirectOverYears"
)
associationByOverallIndirectOverYears_file = (
    results_path + "associationByOverallIndirectOverYears"
)
associationByOverallDirectOverYears_file = (
    results_path + "associationByOverallDirectOverYears"
)
evidenceIndirect_file = results_path + "evidenceIndirectDated"
evidenceDirect_file = results_path + "evidenceDirectDated"


spark = SparkSession.builder.getOrCreate()


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
    diseases = spark.read.parquet(diseases_file).persist()
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


# try different novelty cutoffs to find novel target-disease associations in 2023
for cutoff in range(1, 10):
    cutoff = cutoff / 10
    pairs = (
        spark.read.parquet(associationByDatasourceDirectOverYears_file)
        .filter(
            (F.col("novelty") > 0)
            & (F.col("novelty") >= cutoff)
            & (F.col("year") == 2023)
        )
        .join(
            getTherapeuticAreaForDisease().select(
                "diseaseId", "therapeuticAreaName", "therapeuticArea"
            ),
            "diseaseId",
            "left",
        )
        # filter only indications
        .filter(
            ~F.col("therapeuticAreaName").isin(
                [
                    "measurement",
                    "biological process",
                    "medical procedure",
                    "phenotype",
                    "injury, poisoning or other complication",
                    "animal disease",
                ]
            )
        )
        # filter very generic indications
        .filter(F.col("therapeuticArea") != F.col("diseaseId"))
        .filter(
            ~F.col("diseaseId").isin(["EFO_0000508", "EFO_0010642", "MONDO_0004992"])
        )
        # add novelty data type
        .join(
            spark.createDataFrame(
                getDatatypeForDatasource(),
                schema=["datasourceId", "datatypeId"],
            ),
            "datasourceId",
            "left",
        )
        .select(
            "targetId",
            "diseaseId",
            "novelty",
            "datatypeId",
            "datasourceId",
            "therapeuticAreaName",
        )
        # add disease info
        .join(
            spark.read.parquet(diseases_file).select(
                F.col("id").alias("diseaseId"), F.col("name").alias("diseaseName")
            ),
            "diseaseId",
            "left",
        )
        # add target info
        .join(
            spark.read.parquet(targets_file).select(
                F.col("id").alias("targetId"),
                F.col("biotype").alias("targetBiotype"),
                F.col("approvedSymbol").alias("targetSymbol"),
            ),
            "targetId",
            "left",
        )
        .distinct()
    )

    # counts
    print(cutoff)
    print("novel targets: {:,d}".format(pairs.select("targetId").distinct().count()))
    print(
        "novel target-disease: {:,d}".format(
            pairs.select("targetId", "diseaseId").distinct().count()
        )
    )

    (
        pairs.groupby("datasourceId")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .write.parquet(
            "gs://ot-team/cfalaguera/novels/targetsByDatasource{}".format(cutoff)
        )
    )

    (
        pairs.groupby("datatypeId")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .write.parquet(
            "gs://ot-team/cfalaguera/novels/targetsByDatatype{}".format(cutoff)
        )
    )

    (
        pairs.groupby("targetBiotype")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .write.parquet(
            "gs://ot-team/cfalaguera/novels/targetsByBiotype{}".format(cutoff)
        )
    )

    (
        pairs.groupby("therapeuticAreaName")
        .agg(F.size(F.collect_set("targetId")).alias("targets"))
        .write.parquet(
            "gs://ot-team/cfalaguera/novels/targetsByTherapeuticAreaName{}".format(
                cutoff
            )
        )
    )

    # list of novel targets in the context of diseases
    pairs.groupby("therapeuticAreaName", "targetSymbol").agg(
        F.size(F.collect_set("diseaseId")).alias("nDisease"),
        F.mean("novelty").alias("meanNovelty"),
        F.max("novelty").alias("maxNovelty"),
        F.concat_ws(",", F.collect_set("targetId")).alias("targetIds"),
        F.concat_ws(",", F.collect_set("datasourceId")).alias("datasources"),
        F.concat_ws(";", F.collect_set("diseaseName")).alias("diseaseNames"),
        F.concat_ws(";", F.collect_set("diseaseId")).alias("diseaseIds"),
    ).write.parquet("gs://ot-team/cfalaguera/novels/novels{}".format(cutoff))
