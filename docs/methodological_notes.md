# Methodological notes

## Summary

The available GSE120575 expression matrix is TPM/log-expression rather than raw counts. DESeq2 is therefore not applied directly to the current matrix. The appropriate validation route is patient-level aggregation of B-cell expression, followed by targeted analysis of antigen-presentation genes.

## Purpose

This document explains the main methodological choices and limitations of the melanoma scRNA-seq reanalysis.

The main project finding is an exploratory responder-associated B-cell/APC HLA class II signal. This signal is biologically plausible, but it requires careful interpretation because the available expression matrix is TPM-based and the number of B-cell-informative patients is limited.

## Data type constraint

The available expression matrix is TPM/log-expression, not unnormalised raw integer counts.

DESeq2 is designed for RNA-seq count data. It models raw counts and estimates size factors and dispersion from count-based measurements. Because TPM values have already been normalised, they should not be used directly as DESeq2 input.

For the current matrix, a safer approach is to use aggregated patient-level log-expression values with limma-style modelling or targeted non-parametric tests.

If raw counts were available, a true patient-level pseudobulk workflow using DESeq2, edgeR, or limma-voom would be appropriate.

## Unit of replication

The biological replicate in this project is the patient, not the individual cell.

Single-cell data provide cellular resolution, but cells from the same patient are not independent biological replicates. Treating each cell as an independent sample can inflate the apparent sample size, underestimate variance, and increase the risk of false positive findings.

For that reason, response-group comparisons should be performed at the patient level rather than by directly comparing all cells as independent observations.

## Project-specific validation strategy

The B-cell HLA class II signal should be evaluated with the patient as the unit of replication.

The planned validation logic is:

- Subset B cells.
- Aggregate B-cell expression for each patient.
- Preserve responder versus non-responder metadata at the patient level.
- Focus on a targeted antigen-presentation gene set.
- Use a method appropriate for TPM/log-expression data.
- Interpret the result as exploratory and requiring validation.

Only 12 of 19 patients contributed at least 5 B cells in the B-cell-focused analysis: 7 responders and 5 non-responders. This limited sample size is why a targeted antigen-presentation gene set is preferable to a broad genome-wide differential-expression claim.

A reasonable targeted antigen-presentation gene set includes:

- HLA-DRA
- HLA-DRB1
- HLA-DPA1
- HLA-DPB1
- HLA-DQA1
- HLA-DQB1
- CD74
- CIITA
- HLA-DMA
- HLA-DMB
- HLA-DOA
- HLA-DOB
- CD40
- CD80
- CD86

## Interpretation

The current finding should be described as an exploratory patient-level immune-contexture signal.

It should not be described as a validated biomarker, a causal mechanism, a genome-wide differential-expression result, a TLS result, or a completed multi-omics analysis.

## Recommended future validation

Future validation would ideally include:

- an independent melanoma immunotherapy cohort,
- raw-count single-cell data for true pseudobulk differential-expression analysis,
- spatial transcriptomics or multiplex imaging to evaluate tissue organisation,
- matched tumour-intrinsic data to connect immune-cell states with tumour biology.

## Related documents

- [Spatial/TLS future directions](spatial_tls_future_directions.md)
- [Multi-omics tumour-intrinsic future directions](multiomics_tumour_intrinsic_future_directions.md)

## References

1. Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications.* 2021;12:5692.
2. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology.* 2014;15:550.
3. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research.* 2015;43:e47.
4. Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology.* 2014;15:R29.
