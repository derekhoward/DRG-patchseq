# DRG-patchseq

This repository contains the code and data analysis pipelines for our Patch-seq study on CMi-fibres/sleeping nociceptors of the dorsal root ganglion (DRG).

The analysis scripts should be run in numerical order. Scripts generate outputs that are used by subsequent steps.

## Output Structure

```
results/
├── SCT-integration/           # Integration results from script 01
│   ├── drg_integrated.RDS
│   └── figures/
├── pig_spatial_integration/   # Pig Visium results from scripts 02, 02b
├── human_spatial_integration/ # Human Visium results from scripts 03, 03b
├── figures/
data/
├── summarized_experiments/    # SummarizedExperiment objects for MetaNeighbor
└── processed/                 # Processed gene lists and ortholog mappings
```