# Timeseries and novelty assessment for Open Targets Platform
This repository contains scripts and data files for timestamping biomedical evidence available in the Open Targets Platform, and calculating novelty scores for target-disease associations in the Platform.

# Overview
With the aim of harnessing and expanding the capabilities of the Open Targets Platform, we have developed a method for tracking biomedical evidence over time and using it to discover potential therapeutic targets with emerging novelty signals. We start with a comprehensive effort to retrospectively timestamp millions of pieces of evidence describing relationships between human targets and diseases that Open Targets has collected over the years. We then reconstruct a timeline for each target-disease association in the Open Targets Platform and define a metric to quantify the degree of novelty of a target in the context of a disease based on the current state of knowledge. 

# System Requirements

## Hardware requirements
The scripts require only a standard computer with enough RAM to support the in-memory operations.

## Software requirements

### OS Requirements
This package is supported for macOS. The package has been tested on the following systems:

- macOS: Sequoia (15.0.1)

### Python Dependencies
```
numpy # 1.26.4
pandas # 2.1.0
pyspark # 3.3.2
scipy # 1.11.4
```

# Execution Guide

1) Parse timestamps for Open Targets Platform evidence from original resources:
```
python timestampsParser.py
```

2) Add timestamps to evidence files in Open Targets Platform:
```
python parsePublicationYearGC.py
```

3) Generate timeseries data and assess target-disease associations' novelty:
```
python timeseries.py
```
To demo `timeseries.py` use `data/demo_evidence` as input containing the timestamped evidence for asthma and TSLP association.
To run the scripts for Open Targets Platform 23.06 as a whole we used a Google Cloud machine with the following specifications:
```
gcloud dataproc clusters create cf-timeseries
  --region europe-west1 --zone europe-west1-d
  --single-node --master-machine-type n2-highmem-128
  --master-boot-disk-size 500 --image-version 2.0-debian10
  --project open-targets-eu-dev
```

Run times are specified in the scripts.

# Extra Analysis

- Plot timeseries for a given target-disease association:
```
import plotTimeseries
plotTimeseries.plotTargetDisease(targetId="ENSG00000145777", # TSLP
                                 diseaseId="MONDO_0004979", # asthma
                                 showScore=True, showEvidence=True, showNovelty=True,
                                 img="data/demo_timelines.png")) # overlap of the 3 plots shown below
```
![alt text](https://github.com/opentargets/timeseries/blob/main/data/demo_timelines.png?raw=true)

- List novel target-disease associations in the Open Targets Platform in 2023:
```
python novels.py
```

- Analyse temporal trends in novel drug targets discovery since 2000.
```
python noveltyTrends.py
```

- Analyse the co-ocurrance of novelty peaks across resources of the same type.
```
python noveltyCorrelation.py
```

- General analysis of outputs:
```
analysis.ipynb
```

# License
This project is covered under the Apache 2.0 License.
