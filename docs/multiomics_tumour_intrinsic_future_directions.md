# Multi-omics tumour-intrinsic future directions

## Summary

The current project is a single-cell transcriptomic reanalysis and does not include matched tumour genomics, TMB, neoantigen load, HLA/B2M status, bulk RNA-seq, proteomics, or spatial transcriptomics. Multi-omics integration is therefore a future direction, not a completed analysis.

## Purpose

This document describes how the exploratory B-cell/APC HLA class II signal could be connected to tumour-intrinsic biology in future matched cohorts.

The current analysis identifies an immune-cell expression signal. A stronger translational follow-up would test whether this signal aligns with tumour antigenicity, antigen-presentation integrity, IFN-gamma pathway activity, and spatial immune organisation.

## Current project boundary

The current GSE120575 reanalysis includes:

- single-cell transcriptomic data,
- patient response metadata,
- cell-type annotations,
- patient-level immune-cell summaries.

It does not include:

- matched WES,
- tumour mutational burden,
- neoantigen load,
- HLA genotype,
- B2M alteration status,
- bulk RNA-seq,
- proteomics,
- spatial transcriptomics.

The project should therefore not be described as a multi-omics analysis.

## Rationale for future multi-omics integration

Checkpoint blockade response is influenced by multiple biological layers.

On the tumour side, relevant features may include:

- tumour mutational burden,
- neoantigen load,
- HLA genotype,
- HLA-I/B2M integrity,
- IFN-gamma pathway activity.

On the immune microenvironment side, relevant features may include:

- B-cell HLA-II activity,
- APC activity,
- T-cell state,
- myeloid/APC composition,
- TLS-like spatial organisation.

The current project addresses part of the immune microenvironment side. A future matched cohort could test whether the B-cell/APC HLA-II signal is part of a broader tumour-immune response axis.

## Candidate future data blocks

A future study could integrate several patient-level data blocks.

Single-cell immune features:

- B-cell HLA-II score,
- APC HLA-II score,
- B-cell fraction,
- T-cell state fractions,
- myeloid/APC composition.

Tumour genomics:

- TMB,
- neoantigen load,
- HLA genotype,
- HLA loss,
- B2M alteration,
- JAK1/JAK2 alteration.

Bulk transcriptomics:

- IFN-gamma response,
- cytolytic score,
- antigen-presentation signature,
- checkpoint gene expression.

Spatial pathology:

- TLS density,
- CD20-positive B-cell density,
- CD3/CD8 proximity,
- CD21-positive follicular dendritic cells,
- MECA79-positive high endothelial venules.

## MOFA+ as an unsupervised framework

MOFA+ could be used in a future matched multi-modal cohort to identify latent tumour-immune axes across data modalities.

For example, a latent factor might capture a shared antigen-presentation-high or immune-responsive axis across single-cell, genomic, bulk transcriptomic, and spatial features.

This would be exploratory and unsupervised. It would help identify biological structure across modalities without using response status as the main driver.

## DIABLO as a supervised framework

DIABLO could be used in a future matched multi-omics cohort if the goal is to identify a sparse multi-omics signature that separates responders from non-responders.

This would require matched data blocks measured in the same patients, response status as a categorical outcome, careful feature selection, cross-validation, and independent validation.

Because supervised multi-omics models can overfit in small cohorts, any DIABLO-like analysis would need strict validation and cautious interpretation.

## Interpretation

The current B-cell/APC HLA class II signal should be described as an exploratory immune-cell finding.

A future multi-omics extension would test whether this signal aligns with tumour antigenicity, intact tumour antigen presentation, IFN-gamma responsiveness, spatial TLS-like organisation, and clinical response.

This would help determine whether the B-cell/APC signal is part of a broader tumour-immune response axis.

## Related documents

- [Methodological notes](methodological_notes.md)
- [Spatial/TLS future directions](spatial_tls_future_directions.md)

## References

1. Argelaguet R, Arnol D, Bredikhin D, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology.* 2020;21:111.
2. Budczies J, Kazdal D, Allgäuer M, et al. Tumour mutational burden: clinical utility, challenges and emerging improvements. *Nature Reviews Clinical Oncology.* 2024.
3. Paschen A, Melero I, Ribas A. Central role of the antigen-presentation and interferon-gamma pathways in resistance to immune checkpoint blockade. *Annual Review of Cancer Biology.* 2022.
4. Singh A, Shannon CP, Gautier B, et al. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics.* 2019;35:3055-3062.
