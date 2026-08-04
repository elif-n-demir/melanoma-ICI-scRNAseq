# Tumour–immune validation framework for melanoma immunotherapy response

## Purpose

This document presents a literature-guided conceptual validation framework. No experimental validation was performed as part of the current project.

The purpose of this framework is to translate an exploratory patient-level single-cell observation into a set of testable follow-up questions. It outlines how tumour, lymphoid and myeloid compartments could be investigated across complementary in vitro, molecular, soluble-factor, immunophenotyping and spatial read-outs.

## Starting point: patient-level melanoma scRNA-seq reanalysis

The starting point is an exploratory reanalysis of public pre-treatment melanoma immune checkpoint inhibitor single-cell RNA-seq data from 19 patients.

No major immune-cell proportion difference remained significant after correction for multiple comparisons. However, responder samples showed a stronger patient-level B-cell/APC HLA class II and CD74-associated antigen-presentation signal.

The B-cell-informative subset was small and B-cell numbers were imbalanced between patients. The result is therefore interpreted as a hypothesis-generating observation rather than a validated biomarker or causal mechanism.

The primary analysis also showed a higher median myeloid/macrophage fraction in non-responders than in responders. This was a directional observation and did not remain statistically significant after multiple-testing correction. It should therefore not be interpreted as evidence that myeloid abundance predicts treatment resistance.

## Boundary of the current dataset

The current scRNA-seq dataset was not designed to resolve neutrophil heterogeneity. This framework therefore does not claim that neutrophil states were identified in the original analysis.

Instead, it defines a follow-up strategy for connecting patient-level immune composition and antigen-presentation observations with myeloid-oriented experimental and spatial validation.

The current analysis has several important limitations:

- Cells provide cellular resolution, but patients provide biological replication.
- Cells from the same patient should not be treated as independent biological replicates.
- The available expression matrix is TPM-based rather than raw-count-based and is not suitable for direct DESeq2 analysis.
- Dissociated scRNA-seq does not preserve tissue architecture or spatial cell–cell relationships.
- The original cohort does not contain matched spatial, proteomic or experimental validation data.
- The analysis does not establish neutrophil states, tertiary lymphoid structures, causal tumour–immune interactions or a clinically validated biomarker.

## Why myeloid-oriented validation is a logical extension

In the primary analysis, the median myeloid/macrophage fraction was higher in non-responders than in responders (0.102 versus 0.009). This directional difference did not remain significant after correction for multiple comparisons and is therefore not presented as evidence of a resistance-associated myeloid phenotype.

Nevertheless, the combination of a responder-associated B-cell/APC antigen-presentation signal and a higher directional myeloid/macrophage fraction in non-responders raises a testable follow-up question: whether distinct tumour–myeloid–lymphoid interactions contribute to early immune checkpoint inhibitor response or resistance.

A logical extension would be to examine how melanoma cells, neutrophils or other myeloid populations, and T cells influence one another under controlled experimental conditions. Particular attention could be given to IL-17-associated inflammatory signalling, neutrophil recruitment and activation, antigen-presentation capacity, T-cell activation or dysfunction, and melanoma-cell adaptation during early treatment exposure.


## Experimental validation strategy

### 1. Melanoma–neutrophil–T-cell co-culture models

A staged co-culture design could be used to separate the effects of individual cell compartments before testing more complex interactions.

Initial conditions could include:

- Melanoma cells alone
- Melanoma cells with T cells
- Melanoma cells with neutrophils
- Melanoma cells with both neutrophils and T cells

The same conditions could be compared in conventional 2D culture and in 3D melanoma spheroid or organotypic models. A 2D system would provide a relatively controlled starting point, whereas a 3D model could better represent tumour architecture, cell penetration and local gradients of oxygen, nutrients and soluble mediators.

Experimental conditions could include baseline culture, inflammatory stimulation, and immune checkpoint inhibitor-related treatment conditions where biologically appropriate. Time points should be selected in advance to distinguish early signalling events from later changes in cell viability or phenotype.

The design should also distinguish contact-dependent effects from soluble-factor-mediated effects. This could be approached by comparing direct co-culture with transwell or conditioned-medium experiments.

Key experimental considerations would include:

- Independent biological replicates
- Donor-to-donor variability in primary neutrophils and T cells
- Standardised starting cell numbers and cell ratios
- Cell viability at each collection time point
- Appropriate untreated, single-cell-type and vehicle controls
- Predefined primary and secondary read-outs
- Batch-aware experimental design

Potential read-outs could include melanoma-cell viability and proliferation, immune-cell survival, T-cell activation and cytotoxicity, T-cell dysfunction-associated phenotypes, neutrophil activation or phenotypic-state features, and the relative contribution of direct cell contact and secreted factors.

This section describes a proposed experimental design and does not represent experiments performed as part of the current project.

### 2. Targeted gene-expression read-outs

Targeted gene-expression analysis could be used to test predefined biological questions arising from the co-culture experiments. A focused qPCR panel would be preferable to an unfocused list of markers, with targets selected according to the cell model, stimulation condition and primary hypothesis.

Candidate biological groups could include:

- Neutrophil- and myeloid-associated identity, recruitment and activation markers: `S100A8`, `S100A9`, `FCGR3B`, `CSF3R`, `CXCR2`, `CEACAM8`, `MPO` and `ELANE`
- Interferon-response and T-cell recruitment: `STAT1`, `IRF1`, `CXCL9` and `CXCL10`
- T-cell activation and cytotoxicity: `CD3D`, `CD8A`, `IFNG`, `GZMB` and `PRF1`
- Antigen presentation: `CD74`, `CIITA` and relevant HLA class II genes where supported by the assay and cell type
- Melanoma identity and cell state: `MLANA`, `PMEL`, `TYR`, `SOX10` and `S100B`

These genes represent candidate read-outs rather than a final validated assay panel. The final selection should be narrowed according to the experimental model, literature support, expected expression range and technical assay performance.

A robust qPCR design should include:

- Predefined biological comparisons
- Independent biological replicates
- Technical replicates
- Validation of reference genes under the relevant experimental conditions
- Appropriate negative and positive controls
- Predefined quality-control criteria
- Reporting of effect sizes and uncertainty rather than reliance on p-values alone
- Correction for multiple testing when several targets are evaluated

Gene-expression changes should be interpreted together with cell viability, protein-level or soluble-factor measurements and functional read-outs. A change in a single transcript would not, by itself, establish a neutrophil state, functional immune interaction or mechanism of treatment resistance.

### 3. Cytokine and chemokine profiling

Soluble-factor profiling could be used to determine whether melanoma–neutrophil–T-cell interactions generate inflammatory, immunostimulatory or immunosuppressive signalling patterns that are not captured by transcript measurements alone.

Candidate analytes could include:

- IL-17-axis read-outs: IL-17A and IL-17F, together with related downstream inflammatory mediators where biologically justified
- Neutrophil recruitment and activation: CXCL8/IL-8, G-CSF/CSF3 and GM-CSF/CSF2
- T-cell recruitment and interferon-associated signalling: CXCL9, CXCL10 and IFN-γ
- Myeloid and tumour-associated inflammation: IL-6, TNF, CCL2, CCL3 and CCL4
- Immunoregulatory signals where biologically justified, such as TGF-β

The final panel should be selected according to the experimental model and primary hypothesis rather than by measuring the largest possible number of analytes.

Important design considerations would include:

- Predefined supernatant-collection time points
- Standardised starting cell numbers and culture volumes
- Measurement of cell viability and cell recovery at each time point
- Normalisation strategies that account for differences in viable cell number
- Independent biological replicates and donor variability
- Appropriate medium-only, untreated and single-cell-type controls
- Consistent sample handling, storage and freeze–thaw conditions
- Batch-aware plate design
- Correction for multiple testing when broad multiplex panels are used

Multiplex bead-based assays could provide an efficient initial screen, whereas targeted ELISA or orthogonal protein-level assays could be used to confirm selected findings.

Supportive evidence would require reproducible changes across biological replicates and biologically connected read-outs. An isolated cytokine difference without corresponding cellular, transcriptional or functional evidence would not be sufficient to establish a mechanism of immune checkpoint inhibitor response or resistance.

### 4. Flow cytometry and immunophenotyping

Flow cytometry could be used to quantify how melanoma–neutrophil–T-cell co-culture conditions alter cell survival, activation and functional phenotype.

A staged analysis should first identify the major cellular compartments and then evaluate condition-specific phenotypic changes within each population.

Key technical requirements would include:

- Exclusion of debris and dead cells
- Singlet gating
- Clear separation of melanoma and immune-cell compartments
- Appropriate Fc-receptor blocking
- Compensation or spectral-unmixing controls
- Fluorescence-minus-one controls for markers with uncertain gating boundaries
- Consistent instrument settings and acquisition strategy
- Biological replicates and donor-aware analysis
- Batch controls where experiments are performed on different days

Potential T-cell read-outs could include:

- Viability and recovery
- Activation
- Proliferation
- Cytotoxic potential
- Dysfunction- or exhaustion-associated phenotype

Potential neutrophil or myeloid read-outs could include:

- Identity and purity
- Viability
- Activation
- Degranulation-related phenotype
- Maturation and phenotypic-state features
- Changes induced by direct tumour-cell contact or soluble factors

Melanoma-cell read-outs could include viability, proliferation, immune-mediated killing and treatment-associated phenotypic changes.

Marker panels should be selected according to the experimental question, sample source and instrument configuration. No single surface marker should be treated as definitive proof of a neutrophil state or functional mechanism.

Flow-cytometric findings should be interpreted together with gene-expression, soluble-factor and functional measurements.

### 5. Multiplex immunofluorescence and spatial characterisation

Multiplex imaging could be used to determine not only which tumour and immune-cell populations are present, but also how they are spatially organised within the melanoma microenvironment.

This would address a major limitation of dissociated single-cell RNA-seq, which loses information about tissue architecture, cellular neighbourhoods and physical proximity between melanoma cells, neutrophils, other myeloid populations and T cells.

Potential spatial questions could include:

- Whether neutrophil- or myeloid-rich regions are located within the tumour core, at the invasive margin or in surrounding stromal areas
- Whether neutrophil-related signals are associated with T-cell exclusion, proximity or altered functional phenotype
- Whether antigen-presentation-rich regions overlap with lymphoid or myeloid neighbourhoods
- Whether early treatment exposure changes the spatial organisation of tumour and immune compartments
- Whether distinct spatial patterns are reproducible across independent tissue samples

Candidate marker classes could include:

- Melanoma identity: SOX10, S100B, Melan-A/MLANA and PMEL
- T cells: CD3, CD4, CD8 and FOXP3
- B cells: CD20/MS4A1 and CD79A
- Myeloid and neutrophil-related markers: CD68, CD163, MPO, S100A8, S100A9 and neutrophil elastase
- Antigen presentation: HLA-DR and CD74 where technically supported
- Functional markers: PD-1, PD-L1, Granzyme B and Ki-67

Marker combinations should be selected according to antibody performance, tissue preservation, panel compatibility and the biological question. A single marker should not be treated as definitive proof of cell identity or functional state.

Image analysis could include:

- Cell segmentation and phenotype assignment using predefined marker combinations
- Cell-density measurements within defined tissue compartments
- Nearest-neighbour and distance analyses
- Cell-neighbourhood or spatial-community analysis
- Comparison of tumour-core, invasive-margin and stromal regions where histological annotation is available
- Sample-level statistical analysis rather than treating individual cells as independent biological replicates

Multiplex platforms such as PhenoCycler could support high-dimensional spatial phenotyping, but the resulting cell assignments and neighbourhoods would still require quality control, appropriate pathology review and independent validation.

Spatial proximity alone would not establish direct interaction or causality. Observed neighbourhood patterns should therefore be interpreted together with co-culture experiments, cytokine measurements, gene-expression read-outs and functional assays.

## Computational integration strategy

The experimental and spatial read-outs should be integrated at the level of the biological sample or donor rather than by treating individual cells, image objects or technical measurements as independent biological replicates.

A practical integration strategy could include:

- Harmonised metadata across culture, qPCR, cytokine, flow-cytometry and imaging experiments
- Predefined sample identifiers, treatment conditions and collection time points
- Sample-level summaries for each assay
- Effect sizes and confidence intervals alongside statistical significance
- Explicit modelling of donor and batch effects
- Cross-assay comparison of biologically related features
- Visualisation of concordant and discordant signals across assays
- Independent validation in additional samples where available

For example, a co-culture condition associated with reduced T-cell activation would be more convincing if it were accompanied by compatible cytokine changes, reproducible gene-expression shifts and a spatial pattern suggesting T-cell exclusion or altered tumour–immune neighbourhoods.

Exploratory dimensionality-reduction or clustering methods could be used to visualise integrated patterns, but they should not replace predefined biological comparisons or independent validation.

Advanced multi-omics frameworks such as MOFA+ or DIABLO would only be appropriate if matched data blocks were available from the same biological samples. They were not applied in the current project and should be regarded only as possible future approaches.

The aim of integration would not be to create a clinical predictor from a small exploratory dataset. It would be to identify reproducible tumour–immune relationships that can be tested across complementary experimental systems.

## What would count as supportive evidence

A proposed tumour–immune mechanism would be considered more credible if the same biological direction were supported across several independent read-outs rather than by a single marker or assay.

Examples of supportive evidence could include:

- Reproducible effects across independent biological replicates or donors
- Consistent changes in both 2D and 3D experimental systems
- A co-culture-specific effect that is absent from the corresponding single-cell-type controls
- Concordant gene-expression and protein-level changes
- Cytokine or chemokine changes that are compatible with the observed cellular phenotype
- Reproducible changes in T-cell activation, cytotoxicity or dysfunction
- Neutrophil- or myeloid-associated changes supported by multiple markers and functional read-outs
- Evidence distinguishing direct cell-contact effects from soluble-factor-mediated effects
- Spatial localisation or cell-neighbourhood patterns that are consistent with the experimental observations
- Similar findings in an independent sample set or through an orthogonal assay

For example, a proposed suppressive neutrophil–tumour interaction would be more strongly supported if a defined co-culture condition reproducibly reduced T-cell activation, produced compatible cytokine changes, altered tumour or immune-cell gene expression and was associated with a corresponding spatial pattern in tissue.

Evidence should be evaluated at the sample or donor level, with effect sizes, uncertainty and consistency across replicates considered alongside statistical significance.

No individual read-out would be expected to establish a complete causal mechanism. The strongest support would come from convergence across experimental, molecular, soluble-factor, functional and spatial measurements.

## What would not be sufficient

The following observations would not, on their own, be sufficient to establish neutrophil heterogeneity, a causal tumour–immune mechanism or a clinically relevant marker of immune checkpoint inhibitor response:

- A change in a single marker interpreted as proof of a distinct cell state
- A small number of cells treated as independent biological replicates
- A non-significant directional difference presented as a resistance-associated biomarker
- Expression of one or two neutrophil-related genes used to assign a definitive neutrophil phenotype
- Spatial marker expression interpreted as single-cell annotation without appropriate cell segmentation and validation
- Visium spot-level expression interpreted as proof of individual neutrophil states
- Spatial proximity presented as evidence of direct interaction or causality
- A public-data tutorial presented as validation of the original responder/non-responder cohort
- An HLA class II module calculated after silently excluding unavailable HLA genes
- A cytokine difference interpreted without accounting for cell number, viability, time point or batch
- A co-culture effect interpreted without the corresponding single-cell-type and treatment controls
- Technical replicates presented as independent biological replication
- An exploratory clustering pattern presented as a predefined biological subtype
- A proposed experimental design described as completed wet-lab work
- Concordance within a single small cohort presented as independent validation

The original single-cell dataset provides cellular resolution but not spatial architecture, experimental perturbation or direct evidence of causality. Likewise, a co-culture system can test specific interactions under controlled conditions but cannot fully reproduce the complexity of a patient tumour microenvironment.

Support for a biological mechanism would therefore require convergence across independent biological replicates, complementary assay types and, where possible, an independent validation cohort.

## Limitations

This framework is literature-guided and conceptual. It does not represent completed experimental validation.

The main limitations are:

- The original scRNA-seq dataset was not designed to resolve neutrophil heterogeneity.
- The cohort was small and included only 19 pre-treatment patients.
- The B-cell-informative subset was limited to 12 patients and B-cell numbers were imbalanced.
- No major immune-cell proportion difference remained significant after multiple-testing correction.
- The higher myeloid/macrophage fraction in non-responders was directional and not statistically significant after correction.
- The available expression matrix is TPM-based rather than raw-count-based.
- The analysis does not contain matched spatial tissue, proteomic, metabolic, genomic or cytokine data.
- No 2D or 3D co-culture, qPCR, flow-cytometry, cytokine-profiling, imaging or in vivo experiments were performed.
- Dissociated scRNA-seq does not preserve tissue architecture or direct cell–cell interactions.
- In vitro systems cannot fully reproduce the complexity, treatment history and cellular diversity of a patient tumour microenvironment.
- Primary neutrophils are short-lived and can change phenotype during isolation and culture.
- Donor variability, batch effects and differences in cell viability could influence experimental read-outs.
- Spatial proximity would not, by itself, demonstrate functional interaction or causality.
- Any proposed marker panel would require optimisation and validation for the selected tissue, cell source, assay and instrument.

The framework should therefore be interpreted as a structured route from an exploratory computational observation to testable experimental questions, rather than as evidence that the proposed biological mechanism has already been demonstrated.

## Selected references

1. Sade-Feldman M, et al. Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell*. 2018.
   https://pubmed.ncbi.nlm.nih.gov/30388456/

2. Squair JW, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021.
   https://www.nature.com/articles/s41467-021-25960-2

3. Jaillon S, et al. Neutrophil diversity and plasticity in tumour progression and therapy. *Nature Reviews Cancer*. 2020.
   https://www.nature.com/articles/s41568-020-0281-y

4. Hirschhorn D, et al. T cell immunotherapies engage neutrophils to eliminate tumor antigen escape variants. *Cell*. 2023.
   https://pubmed.ncbi.nlm.nih.gov/37001503/

5. Schlenker R, et al. Myeloid–T cell interplay and cell state transitions associated with checkpoint inhibitor response in melanoma. *Med*. 2024.
   https://pubmed.ncbi.nlm.nih.gov/38593812/

6. Cabrita R, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature*. 2020.
   https://www.nature.com/articles/s41586-019-1914-8

7. Helmink BA, et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature*. 2020.
   https://www.nature.com/articles/s41586-019-1922-8

8. Di Blasio S, et al. The tumour microenvironment shapes dendritic cell plasticity in a human organotypic melanoma culture. *Nature Communications*. 2020.
   https://www.nature.com/articles/s41467-020-16583-0

9. Antoranz A, Van Herck Y, Bolognesi MM, et al. Mapping the Immune Landscape in Metastatic Melanoma Reveals Localized Cell–Cell Interactions That Predict Immunotherapy Response. *Cancer Research*. 2022;82(18):3275–3290.
   https://pubmed.ncbi.nlm.nih.gov/35834277/