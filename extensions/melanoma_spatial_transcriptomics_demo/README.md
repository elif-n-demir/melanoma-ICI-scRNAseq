# Melanoma spatial transcriptomics demo

This extension is an independent workflow demonstration using a public human melanoma spatial transcriptomics dataset. It is not a spatial validation of the GSE120575 responder/non-responder single-cell RNA-seq analysis presented in the main repository.

## Purpose

The purpose of this demo is to show a basic Seurat workflow for loading, processing and visualising spatial gene-expression data from an independent human melanoma sample.

The analysis focuses on selected melanoma-associated, immune-cell-associated, myeloid-related and antigen-presentation-related markers. It is intended as a technical and exploratory workflow demonstration, not as a biological validation study.

## Dataset

This demo uses a public 10x Genomics human melanoma FFPE spatial gene-expression dataset generated with the Visium CytAssist platform.

The raw data are stored locally under:

`data_raw/10x_human_melanoma_if_ffpe/`

Raw data are not committed to GitHub.

The required local input files are:

- `CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5`
- `spatial/`

## Workflow

The analysis script performs the following steps:

1. Loads the spatial dataset with `Load10X_Spatial()`.
2. Calculates basic spot-level quality-control metrics.
3. Applies standard log-normalisation.
4. Identifies variable features.
5. Scales the data and runs PCA.
6. Performs neighbour detection, clustering and UMAP.
7. Generates a spatial cluster plot.
8. Visualises selected marker genes on the tissue section.
9. Calculates a targeted antigen-presentation module score.
10. Saves figures, result tables and session information.

The analysis script is located at:

`scripts/01_melanoma_visium_spatial_demo.R`

## Selected marker groups

Melanoma-associated markers:

- SOX10
- MLANA
- PMEL
- S100B

T-cell-associated markers:

- CD3D
- CD3E
- CD8A

Antigen-presentation-associated markers:

- CD74
- HLA-DRA

Myeloid/neutrophil-related markers:

- S100A8
- S100A9

These markers are used for exploratory spatial visualisation only. Marker expression in Visium data should not be interpreted as definitive single-cell annotation.

## Outputs

Figures are saved in:

`figures/`

Result tables are saved in:

`results/`

Main figure outputs include:

- `01_qc_violin_plots.png`
- `02_spatial_clusters.png`
- `03_spatial_marker_featureplots.png`
- `04_antigen_presentation_module_score.png`

Main result outputs include:

- `01_object_summary.csv`
- `03_marker_gene_presence.csv`
- `04_antigen_presentation_gene_presence.csv`
- `05_antigen_presentation_module_score_by_spot.csv`
- `07_sessionInfo.txt`

## Important limitations

This demo has several important limitations.

It uses an independent spatial transcriptomics dataset, not the GSE120575 melanoma immune-checkpoint-inhibitor cohort. It does not contain responder/non-responder metadata and does not validate the responder-associated B-cell/APC HLA class II signal identified in the main single-cell analysis.

Visium data are spot-level data and should not be interpreted as single-cell-resolution data. Myeloid- or neutrophil-related marker expression does not prove neutrophil heterogeneity. Spatial proximity does not prove direct cell-cell interaction or causality.

The antigen-presentation module score is exploratory and should not be interpreted as a validated biomarker.

## Relationship to the main repository

The main repository presents an exploratory reanalysis of public pre-treatment melanoma immune checkpoint inhibitor single-cell RNA-seq data.

This spatial transcriptomics demo is included only as an additional technical extension to show how an independent human melanoma spatial dataset can be loaded, processed and visualised in Seurat.

It should be interpreted as a workflow demonstration, not as a validation experiment.