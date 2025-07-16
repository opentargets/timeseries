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
requests # 2.31.0
```

# Execution Guide

1) Parse timestamps for Open Targets Platform evidence from original resources:
```
python parse_timestamp.py # run locally (ERT: 3 min) and output files are automatically uploaded to GC
```

2) Setup Google Cloud machine:
```
gcloud dataproc clusters create cf-timeseries
  --image-version 2.2 --region europe-west1
  --master-machine-type n1-standard-2
  --secondary-worker-type spot --worker-machine-type n1-standard-4 --worker-boot-disk-size 500
  --autoscaling-policy=otg-etl --enable-component-gateway
  --project open-targets-eu-dev
```

3) Add timestamps to evidence files in Open Targets Platform:
```
python timestamp_evidence.py # run in GC (ERT: 6 min)
```

3) Generate timeseries data and assess target-disease associations' novelty:
```
python timeseries.py # run in GC (ERT: 30 min)
```

# Extra Analysis

- Plot timeseries for a given target-disease association:
```
# run locally
import plot_timeseries
plot_timeseries.plotTargetDisease(targetId="ENSG00000145777", # TSLP
                                 diseaseId="MONDO_0004979", # asthma
                                 showScore=True, showEvidence=True, showNovelty=True,
                                 img="data/demo_timelines.png")) # overlap of the 3 plots shown below
```
![alt text](https://github.com/opentargets/timeseries/blob/main/data/demo_timelines.png?raw=true)

- Count the number of novel target-disease associations and unique novel targets over the years across Open Targets Platform resources.
```
python novelty_distribution.py # run in GC
```

- Analyse temporal trends in novel drug targets discovery from 2000 to 2025.
```
python novelty_approval.py # run in GC
```

- General analysis of outputs:
```
analysis.ipynb
```

# Data availability

The entire data generated with this code has been deposited on Zenodo: https://zenodo.org/records/15922783.

# License
This project is covered under the Apache 2.0 License.
