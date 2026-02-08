# DRG-patchseq

This repository contains the code and data analysis pipelines for our Patch-seq study on CMi-fibres/sleeping nociceptors of the dorsal root ganglion (DRG).

Körner, J.\*, Howard, D.\*, et al., Tripathy, S.† & Lampert, A.† (2026).  
*Molecular architecture of human dermal sleeping nociceptors*. Cell.  
\*Co–first authors. †Corresponding authors.

## Data availability
Gene expression data is available from [GEO:GSE263532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263532)

The data from the harmonized cross-species atlas used in this study is also available from [GEO:GSE255436](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255436)

## Output Structure
The analysis scripts should be run in numerical order. Scripts generate outputs that are used by subsequent steps.

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