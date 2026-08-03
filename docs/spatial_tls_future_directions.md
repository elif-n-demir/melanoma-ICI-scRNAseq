# Spatial/TLS future directions

## Summary

The B-cell/APC HLA class II signal is compatible with B-cell/TLS-related immunotherapy biology, but single-cell RNA-seq alone cannot prove tertiary lymphoid structures. Spatial or histological validation would be required before making any TLS claim.

## Purpose

This document describes how the exploratory B-cell/APC HLA class II signal could be validated spatially in future work.

The current project is based on dissociated single-cell RNA-seq. This approach can identify cell states and gene-expression patterns, but it does not preserve tissue architecture. It cannot determine whether B cells, T cells, follicular dendritic cells, and high endothelial venules are spatially organised into tertiary lymphoid structures.

## Current project boundary

The current analysis supports an exploratory responder-associated B-cell/APC HLA class II antigen-presentation signal.

The current analysis does not identify tertiary lymphoid structures.

TLS require spatial or histological evidence. The scRNA-seq result should therefore be interpreted as hypothesis-generating, not as direct TLS evidence.

## Spatial validation objective

A future spatial validation study would assess whether HLA-II-high B cells are organised near T cells within TLS-like immune niches.

The main objective would be to evaluate whether HLA-II-high B-cell regions in responder samples overlap with T-cell-rich areas and TLS-associated structures.

## Candidate validation methods

### Multiplex immunofluorescence or multiplex IHC

A multiplex staining panel could include several marker groups.

B-cell markers:

- CD20
- CD19
- MS4A1

HLA-II and antigen-presentation markers:

- HLA-DR
- HLA-DRA
- CD74

T-cell markers:

- CD3
- CD4
- CD8
- FOXP3

TLS-associated markers:

- CXCL13
- CXCR5
- CD21
- CD23
- MECA79
- LAMP3 / DC-LAMP

This approach would allow assessment of whether B cells, T cells, and TLS-associated structures are co-localised in the tissue.

### Spatial transcriptomics

Spatial transcriptomics could measure gene expression while preserving spatial coordinates within the tissue section. This would allow evaluation of whether HLA-II-high B-cell regions overlap with T-cell-rich or TLS-like areas.

Public melanoma spatial transcriptomic data have been reported, including the Thrane et al. study of stage III cutaneous melanoma lymph node metastases. The number of melanoma spatial datasets is also increasing, making this a realistic future validation direction.

### Digital spatial profiling

Region-of-interest-based spatial profiling could compare immune markers inside and outside TLS-like regions. This could help evaluate whether T cells near B-cell-rich regions show different activation or exhaustion patterns.

### Pathology review

Histological review using H&E and IHC would be needed to distinguish random lymphocyte aggregates from organised TLS-like structures.

## Biological interpretation

The B-cell/APC HLA class II signal is consistent with literature linking B cells and TLS-related immune organisation to immunotherapy response in melanoma. However, it remains an expression-based signal. It should be interpreted as a hypothesis-generating finding rather than direct evidence of TLS.

## Recommended future workflow

A future spatial validation workflow could include:

- selecting pre-treatment melanoma samples with known response status,
- staining for B-cell, T-cell, HLA-II, and TLS-associated markers,
- quantifying spatial proximity between B cells and T cells,
- assessing TLS morphology and maturity,
- comparing spatial features between responders and non-responders,
- integrating spatial findings with the single-cell HLA-II signal.

## Related documents

- [Methodological notes](methodological_notes.md)
- [Multi-omics tumour-intrinsic future directions](multiomics_tumour_intrinsic_future_directions.md)

## References

1. Cabrita R, Lauss M, Sanna A, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature.* 2020;577:561-565.
2. Helmink BA, Reddy SM, Gao J, et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature.* 2020;577:549-555.
3. Petitprez F, de Reyniès A, Keung EZ, et al. B cells are associated with survival and immunotherapy response in sarcoma. *Nature.* 2020;577:556-560.
4. Thrane K, Eriksson H, Maaskola J, Hansson J, Lundeberg J. Spatially resolved transcriptomics enables dissection of genetic heterogeneity in stage III cutaneous malignant melanoma. *Cancer Research.* 2018;78:5970-5979.
