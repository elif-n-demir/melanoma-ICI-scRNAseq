# =============================================================================
# 06_melanoma_qc_normalization.R
# Purpose: QC and normalization for GSE120575 pre-treatment melanoma object
# =============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("data_processed", showWarnings = FALSE)

mel <- readRDS("data_processed/GSE120575_pre_seurat_raw.rds")

mel
table(mel$response)
table(mel$patient_id)

mel[["percent.mt"]] <- PercentageFeatureSet(
  mel,
  pattern = "^MT-"
)

summary(mel$nFeature_RNA)
summary(mel$nCount_RNA)
summary(mel$percent.mt)

p_qc_patient <- VlnPlot(
  mel,
  features = c("nFeature_RNA", "percent.mt"),
  group.by = "patient_id",
  pt.size = 0,
  ncol = 2
)

p_qc_patient

ggsave(
  filename = "figures/GSE120575_pre_QC_by_patient.png",
  plot = p_qc_patient,
  width = 14,
  height = 6,
  dpi = 300
)

# =============================================================================
# Normalization and variable features
# =============================================================================

mel <- NormalizeData(
  mel,
  normalization.method = "LogNormalize"
)

mel <- FindVariableFeatures(
  mel,
  selection.method = "vst",
  nfeatures = 2000
)

top20_variable_genes <- head(
  VariableFeatures(mel),
  20
)

top20_variable_genes

p_variable_features <- VariableFeaturePlot(mel)

p_variable_features_labeled <- LabelPoints(
  plot = p_variable_features,
  points = top20_variable_genes,
  repel = TRUE
)

p_variable_features_labeled

ggsave(
  filename = "figures/GSE120575_pre_variable_features.png",
  plot = p_variable_features_labeled,
  width = 8,
  height = 6,
  dpi = 300
)

write_csv(
  tibble(gene = top20_variable_genes),
  "results/GSE120575_pre_top20_variable_genes.csv"
)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_normalized_variablefeatures.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_normalized_variablefeatures.rds")

# =============================================================================
# Scaling and PCA
# =============================================================================

mel <- ScaleData(
  mel,
  features = VariableFeatures(mel)
)

mel <- RunPCA(
  mel,
  features = VariableFeatures(mel)
)

mel[["pca"]]

p_pca <- DimPlot(
  mel,
  reduction = "pca",
  group.by = "response"
)

p_pca

ggsave(
  filename = "figures/GSE120575_pre_PCA_by_response.png",
  plot = p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

p_elbow <- ElbowPlot(
  mel,
  ndims = 50
)

p_elbow

ggsave(
  filename = "figures/GSE120575_pre_elbow_plot.png",
  plot = p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

# =============================================================================
# Clustering and UMAP
# =============================================================================

dims_use <- 1:20

mel <- FindNeighbors(
  mel,
  dims = dims_use
)

mel <- FindClusters(
  mel,
  resolution = 0.4
)

mel <- RunUMAP(
  mel,
  dims = dims_use
)

mel

table(mel$seurat_clusters)

p_umap_clusters <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

p_umap_response <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "response"
)

p_umap_patient <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "patient_id"
)

p_umap_clusters
p_umap_response
p_umap_patient

# =============================================================================
# Save UMAP plots and clustered object
# =============================================================================

ggsave(
  filename = "figures/GSE120575_pre_UMAP_clusters_res0.4.png",
  plot = p_umap_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_UMAP_response.png",
  plot = p_umap_response,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_UMAP_patient.png",
  plot = p_umap_patient,
  width = 8,
  height = 6,
  dpi = 300
)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_clustered_res0.4.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_clustered_res0.4.rds")

# =============================================================================
# Cluster marker analysis
# =============================================================================

mel_markers <- FindAllMarkers(
  mel,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

dim(mel_markers)

head(mel_markers)

top10_markers <- mel_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

top10_markers %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = 200)

write_csv(
  mel_markers,
  "results/GSE120575_pre_all_cluster_markers.csv"
)

write_csv(
  top10_markers,
  "results/GSE120575_pre_top10_cluster_markers.csv"
)

file.exists("results/GSE120575_pre_all_cluster_markers.csv")
file.exists("results/GSE120575_pre_top10_cluster_markers.csv")

melanoma_marker_panel <- c("MLANA", "PMEL", "MITF", "TYR")

p_melanoma_control <- DotPlot(
  mel,
  features = melanoma_marker_panel
) +
  RotatedAxis()

p_melanoma_control

melanoma_control_table <- p_melanoma_control$data %>%
  select(
    cluster = id,
    gene = features.plot,
    avg.exp,
    pct.exp,
    avg.exp.scaled
  ) %>%
  arrange(
    gene,
    desc(pct.exp)
  )

melanoma_control_table <- as_tibble(melanoma_control_table)

print(
  melanoma_control_table,
  n = 40
)

# =============================================================================
# Check whether MITF-high cluster 6 is myeloid/macrophage-like or melanoma-like
# =============================================================================

myeloid_melanoma_check_panel <- c(
  "MLANA", "PMEL", "MITF", "TYR",
  "MARCO", "OLR1", "CD300E", "MAFB", "FPR3",
  "LYZ", "CD14", "FCGR3A"
)

myeloid_melanoma_check_panel_present <- myeloid_melanoma_check_panel[
  myeloid_melanoma_check_panel %in% rownames(mel)
]

myeloid_melanoma_check_panel_present

p_myeloid_melanoma_check <- DotPlot(
  mel,
  features = myeloid_melanoma_check_panel_present
) +
  RotatedAxis()

p_myeloid_melanoma_check

melanoma_control_table %>%
  print(n = 100)

# =============================================================================
# Save melanoma / myeloid control plots
# =============================================================================

ggsave(
  filename = "figures/GSE120575_pre_melanoma_marker_control_dotplot.png",
  plot = p_melanoma_control,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_myeloid_melanoma_check_dotplot.png",
  plot = p_myeloid_melanoma_check,
  width = 10,
  height = 5,
  dpi = 300
)

write_csv(
  melanoma_control_table,
  "results/GSE120575_pre_melanoma_marker_control_table.csv"
)

# =============================================================================
# Immune cell type annotation marker panel
# =============================================================================

annotation_marker_panel <- c(
  # Pan T cells
  "CD3D", "CD3E", "TRAC",
  
  # CD4 / naive-memory T
  "CD4", "IL7R", "CCR7", "TCF7", "LEF1",
  
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "ICOS",
  
  # CD8 / cytotoxic T
  "CD8A", "CD8B", "GZMK", "GZMB", "PRF1",
  
  # NK / cytotoxic
  "NKG7", "GNLY", "KLRD1", "KLRF1", "FGFBP2",
  
  # B cells
  "MS4A1", "CD79A", "CD79B", "CD22", "BANK1",
  
  # Plasma cells
  "SDC1", "MZB1", "JCHAIN", "IGKC",
  
  # Myeloid / macrophage
  "LYZ", "CD14", "FCGR3A", "MARCO", "OLR1", "MAFB",
  
  # Dendritic / pDC
  "FCER1A", "CLEC10A", "CLEC4C", "LILRA4",
  
  # Mast cell
  "TPSAB1", "CPA3"
)

annotation_marker_panel_present <- annotation_marker_panel[
  annotation_marker_panel %in% rownames(mel)
]

annotation_marker_panel_present

p_annotation_markers <- DotPlot(
  mel,
  features = annotation_marker_panel_present
) +
  RotatedAxis()

p_annotation_markers

ggsave(
  filename = "figures/GSE120575_pre_annotation_marker_dotplot.png",
  plot = p_annotation_markers,
  width = 16,
  height = 7,
  dpi = 300
)

# =============================================================================
# First-pass cell type annotation
# =============================================================================

cluster_annotations <- c(
  "0" = "Treg / activated CD4 T",
  "1" = "Activated cytotoxic CD8 T",
  "2" = "Naive-memory T",
  "3" = "Cycling activated T",
  "4" = "NK / cytotoxic T",
  "5" = "B cell",
  "6" = "Myeloid / macrophage",
  "7" = "Plasma cell",
  "8" = "Ig-high B / plasma-like",
  "9" = "pDC / dendritic-like"
)

mel$cell_type <- cluster_annotations[
  as.character(mel$seurat_clusters)
]

table(mel$seurat_clusters, mel$cell_type)

cluster_annotations <- c(
  "0" = "Treg / activated CD4 T",
  "1" = "Activated cytotoxic CD8 T",
  "2" = "Naive-memory T",
  "3" = "Cycling activated T",
  "4" = "NK / cytotoxic T",
  "5" = "B cell",
  "6" = "Myeloid / macrophage",
  "7" = "Plasma cell",
  "8" = "Ig-high B / plasma-like",
  "9" = "pDC / dendritic-like"
)

cell_type_vector <- cluster_annotations[
  as.character(mel$seurat_clusters)
]

cell_type_df <- data.frame(
  cell_type = unname(cell_type_vector),
  row.names = colnames(mel)
)

mel <- AddMetaData(
  object = mel,
  metadata = cell_type_df
)

table(mel$seurat_clusters, mel$cell_type)

p_umap_celltype <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)

p_umap_celltype

# =============================================================================
# Save annotated UMAP and annotated Seurat object
# =============================================================================

ggsave(
  filename = "figures/GSE120575_pre_UMAP_celltype_annotation.png",
  plot = p_umap_celltype,
  width = 10,
  height = 8,
  dpi = 300
)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_annotated.rds"
)

file.exists("figures/GSE120575_pre_UMAP_celltype_annotation.png")
file.exists("data_processed/GSE120575_pre_seurat_annotated.rds")

celltype_counts <- mel@meta.data %>%
  count(cell_type, response, name = "n_cells") %>%
  arrange(cell_type, response)

celltype_counts

write_csv(
  celltype_counts,
  "results/GSE120575_pre_celltype_counts_by_response.csv"
)

# =============================================================================
# T cell annotation refinement: CD4 / CD8 / Treg panel
# =============================================================================

tcell_cd4_cd8_panel <- c(
  "CD3D", "CD3E", "TRAC",
  "CD4", "CD8A", "CD8B",
  "IL7R", "CCR7", "TCF7", "LEF1",
  "FOXP3", "IL2RA", "CTLA4", "ICOS"
)

tcell_cd4_cd8_panel_present <- tcell_cd4_cd8_panel[
  tcell_cd4_cd8_panel %in% rownames(mel)
]

tcell_cd4_cd8_panel_present

p_tcell_cd4_cd8_dotplot <- DotPlot(
  mel,
  features = tcell_cd4_cd8_panel_present
) +
  RotatedAxis()

p_tcell_cd4_cd8_dotplot

p_tcell_cd4_cd8_vln <- VlnPlot(
  mel,
  features = c("CD4", "CD8A", "CD8B", "IL7R", "FOXP3"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)

p_tcell_cd4_cd8_vln

# =============================================================================
# T cell exhaustion / tumor-reactive marker panel
# =============================================================================

exhaustion_panel <- c(
  "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1",
  "CTLA4", "CXCL13", "TNFRSF9", "GZMK",
  "TCF7", "SELL", "CCR7", "IL7R"
)

exhaustion_panel_present <- exhaustion_panel[
  exhaustion_panel %in% rownames(mel)
]

exhaustion_panel_present

p_exhaustion_dotplot <- DotPlot(
  mel,
  features = exhaustion_panel_present
) +
  RotatedAxis()

p_exhaustion_dotplot

# =============================================================================
# Refined cell type annotation after T cell marker checks
# =============================================================================

cluster_annotations_refined <- c(
  "0" = "Treg",
  "1" = "GZMK+ activated/transitional CD8 T",
  "2" = "TCF7+ memory-like CD8 T",
  "3" = "Cycling exhausted/tumor-reactive CD8 T",
  "4" = "NK / cytotoxic T",
  "5" = "B cell",
  "6" = "Myeloid / macrophage",
  "7" = "Plasma cell",
  "8" = "Plasma cell / Ig-high",
  "9" = "pDC / dendritic-like"
)

cell_type_refined_vector <- cluster_annotations_refined[
  as.character(mel$seurat_clusters)
]

cell_type_refined_df <- data.frame(
  cell_type_refined = unname(cell_type_refined_vector),
  row.names = colnames(mel)
)

mel <- AddMetaData(
  object = mel,
  metadata = cell_type_refined_df
)

table(mel$seurat_clusters, mel$cell_type_refined)

p_umap_celltype_refined <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "cell_type_refined",
  label = TRUE,
  repel = TRUE
)

p_umap_celltype_refined

ggsave(
  filename = "figures/GSE120575_pre_tcell_CD4_CD8_Treg_dotplot.png",
  plot = p_tcell_cd4_cd8_dotplot,
  width = 11,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_tcell_CD4_CD8_Treg_violin.png",
  plot = p_tcell_cd4_cd8_vln,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_tcell_exhaustion_dotplot.png",
  plot = p_exhaustion_dotplot,
  width = 12,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_UMAP_celltype_refined_annotation.png",
  plot = p_umap_celltype_refined,
  width = 10,
  height = 8,
  dpi = 300
)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_refined_annotated.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_refined_annotated.rds")

# =============================================================================
# CD4 / Treg check before proportion analysis
# =============================================================================

Idents(mel) <- "seurat_clusters"

cd4_treg_genes <- c(
  "CD3D", "CD3E", "TRAC",
  "CD4", "CD8A", "CD8B",
  "IL7R", "CCR7", "TCF7", "LEF1",
  "FOXP3", "IL2RA", "CTLA4", "ICOS",
  "CD40LG"
)

cd4_treg_genes_present <- cd4_treg_genes[
  cd4_treg_genes %in% rownames(mel)
]

cd4_treg_genes_present

p_cd4_treg_check <- DotPlot(
  mel,
  features = cd4_treg_genes_present
) +
  RotatedAxis()

p_cd4_treg_check

c0_cells <- WhichCells(
  mel,
  idents = "0"
)

data_mat <- GetAssayData(
  mel,
  assay = "RNA",
  layer = "data"
)

cluster0_marker_positive <- tibble(
  gene = c("FOXP3", "IL2RA", "CTLA4", "ICOS", "CD40LG", "IL7R", "CD4", "CD8A", "CD8B")
) %>%
  filter(gene %in% rownames(data_mat)) %>%
  mutate(
    pct_positive_cluster0 = map_dbl(
      gene,
      ~ mean(data_mat[.x, c0_cells] > 0) * 100
    )
  )

cluster0_marker_positive

# =============================================================================
# Safe resolution 0.8 test to see whether CD4/Treg cluster separates
# =============================================================================

mel_res08 <- mel

mel_res08$old_cluster_res04 <- as.character(mel_res08$seurat_clusters)

mel_res08 <- FindClusters(
  mel_res08,
  resolution = 0.8
)

table(
  old_cluster = mel_res08$old_cluster_res04,
  new_cluster = mel_res08$seurat_clusters
)

p_res08_cd4_treg_check <- DotPlot(
  mel_res08,
  features = cd4_treg_genes_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

p_res08_cd4_treg_check

# =============================================================================
# NK / gamma-delta T check for cluster 4
# =============================================================================

nk_gd_panel <- c(
  "CD3D", "CD3E", "TRAC",
  "TRDC", "TRGC1", "TRGC2",
  "KLRF1", "KLRD1", "NCAM1",
  "FGFBP2", "GNLY", "NKG7", "FCGR3A"
)

nk_gd_panel_present <- nk_gd_panel[
  nk_gd_panel %in% rownames(mel)
]

nk_gd_panel_present

p_nk_gd_check <- DotPlot(
  mel,
  features = nk_gd_panel_present
) +
  RotatedAxis()

p_nk_gd_check

# =============================================================================
# Final cell type annotation after CD4/Treg and NK/gamma-delta checks
# =============================================================================

cluster_annotations_final <- c(
  "0" = "CD4 T / Treg-enriched",
  "1" = "GZMK+ activated/transitional CD8 T",
  "2" = "TCF7+ memory-like CD8 T",
  "3" = "Cycling exhausted/tumor-reactive CD8 T",
  "4" = "gamma-delta / NK-like cytotoxic T",
  "5" = "B cell",
  "6" = "Myeloid / macrophage",
  "7" = "Plasma cell",
  "8" = "Plasma cell / Ig-high",
  "9" = "pDC / dendritic-like"
)

cell_type_final_vector <- cluster_annotations_final[
  as.character(mel$seurat_clusters)
]

cell_type_final_df <- data.frame(
  cell_type_final = unname(cell_type_final_vector),
  row.names = colnames(mel)
)

mel <- AddMetaData(
  object = mel,
  metadata = cell_type_final_df
)

table(mel$seurat_clusters, mel$cell_type_final)

p_umap_celltype_final <- DimPlot(
  mel,
  reduction = "umap",
  group.by = "cell_type_final",
  label = TRUE,
  repel = TRUE
)

p_umap_celltype_final

ggsave(
  filename = "figures/GSE120575_pre_CD4_Treg_check_dotplot.png",
  plot = p_cd4_treg_check,
  width = 11,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_resolution08_CD4_Treg_check_dotplot.png",
  plot = p_res08_cd4_treg_check,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_NK_gamma_delta_check_dotplot.png",
  plot = p_nk_gd_check,
  width = 11,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "figures/GSE120575_pre_UMAP_celltype_final_annotation.png",
  plot = p_umap_celltype_final,
  width = 10,
  height = 8,
  dpi = 300
)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_final_annotated.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_final_annotated.rds")

# =============================================================================
# Final annotation metadata check
# =============================================================================

table(mel$cell_type_final, useNA = "ifany")

sum(is.na(mel$cell_type_final))

# =============================================================================
# Analysis-level cell type labels
# Plasma cell clusters are merged for patient-level proportion analysis
# =============================================================================

mel$cell_type_analysis <- ifelse(
  mel$cell_type_final %in% c("Plasma cell", "Plasma cell / Ig-high"),
  "Plasma cell",
  mel$cell_type_final
)

table(mel$cell_type_final, mel$cell_type_analysis)

# =============================================================================
# Subcluster cluster 0: can CD4 conventional and Treg-like cells be separated?
# =============================================================================

Idents(mel) <- "seurat_clusters"

names(mel@graphs)

mel <- FindSubCluster(
  mel,
  cluster = "0",
  graph.name = "RNA_snn",
  subcluster.name = "sub0",
  resolution = 0.3
)

table(mel$sub0, useNA = "ifany")

# Cluster 0 only
Idents(mel) <- "seurat_clusters"

c0_cells <- WhichCells(
  mel,
  idents = "0"
)

mel_c0 <- subset(
  mel,
  cells = c0_cells
)

sub0_genes <- c(
  "FOXP3", "IL2RA", "CTLA4", "ICOS",
  "CD40LG", "IL7R", "TCF7", "LEF1",
  "CD4", "CD8A", "CD8B"
)

sub0_genes_present <- sub0_genes[
  sub0_genes %in% rownames(mel_c0)
]

p_sub0_cd4_treg <- DotPlot(
  mel_c0,
  features = sub0_genes_present,
  group.by = "sub0"
) +
  RotatedAxis()

p_sub0_cd4_treg

# Percent positive per subcluster inside old cluster 0

data_mat <- GetAssayData(
  mel,
  assay = "RNA",
  layer = "data"
)

sub0_marker_positive <- map_dfr(
  sort(unique(mel_c0$sub0)),
  function(sc) {
    cells_sc <- colnames(mel_c0)[mel_c0$sub0 == sc]
    
    tibble(
      sub0 = sc,
      gene = sub0_genes_present,
      pct_positive = map_dbl(
        sub0_genes_present,
        ~ mean(data_mat[.x, cells_sc] > 0) * 100
      )
    )
  }
)

sub0_marker_positive

# =============================================================================
# Full numeric check of cluster 0 subclusters
# =============================================================================

data_mat <- GetAssayData(
  mel,
  assay = "RNA",
  layer = "data"
)

sub0_genes <- c(
  "FOXP3", "IL2RA", "CTLA4", "ICOS",
  "CD40LG", "IL7R", "TCF7", "LEF1",
  "CD4", "CD8A", "CD8B"
)

sub0_genes_present <- sub0_genes[
  sub0_genes %in% rownames(data_mat)
]

sub0_levels <- sort(unique(mel$sub0[mel$seurat_clusters == "0"]))

sub0_marker_positive <- map_dfr(
  sub0_levels,
  function(sc) {
    cells_sc <- colnames(mel)[
      mel$seurat_clusters == "0" &
        mel$sub0 == sc
    ]
    
    tibble(
      sub0 = sc,
      n_cells = length(cells_sc),
      gene = sub0_genes_present,
      pct_positive = map_dbl(
        sub0_genes_present,
        ~ mean(data_mat[.x, cells_sc] > 0) * 100
      )
    )
  }
)

sub0_marker_positive_wide <- sub0_marker_positive %>%
  select(sub0, n_cells, gene, pct_positive) %>%
  pivot_wider(
    names_from = gene,
    values_from = pct_positive
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

print(sub0_marker_positive_wide, n = Inf)

write_csv(
  sub0_marker_positive,
  "results/GSE120575_pre_cluster0_subcluster_marker_positive.csv"
)

# =============================================================================
# Create analysis-level cell type labels
# =============================================================================

treg_subclusters <- sub0_marker_positive %>%
  filter(gene == "FOXP3", pct_positive >= 50) %>%
  pull(sub0)

treg_subclusters

mel$cell_type_analysis <- mel$cell_type_final

# Merge plasma cell clusters for analysis
mel$cell_type_analysis[
  mel$cell_type_final %in% c("Plasma cell", "Plasma cell / Ig-high")
] <- "Plasma cell"

# Refine cluster 0 only if a clearly FOXP3-high subcluster exists
if (length(treg_subclusters) > 0) {
  
  mel$cell_type_analysis[
    mel$seurat_clusters == "0" &
      mel$sub0 %in% treg_subclusters
  ] <- "Treg-like CD4 T"
  
  mel$cell_type_analysis[
    mel$seurat_clusters == "0" &
      !(mel$sub0 %in% treg_subclusters)
  ] <- "Conventional / memory-like CD4 T"
  
} else {
  
  mel$cell_type_analysis[
    mel$seurat_clusters == "0"
  ] <- "CD4 T / Treg-enriched"
}

table(mel$cell_type_final, mel$cell_type_analysis, useNA = "ifany")
sum(is.na(mel$cell_type_analysis))

# =============================================================================
# Patient-level cell type proportions: all immune cells
# =============================================================================

patient_info <- mel@meta.data %>%
  distinct(patient_id, response)

all_cell_types <- sort(unique(mel$cell_type_analysis))

patient_cell_grid <- patient_info %>%
  tidyr::crossing(cell_type_analysis = all_cell_types)

celltype_counts_all <- mel@meta.data %>%
  count(patient_id, response, cell_type_analysis, name = "n_cells")

prop_all <- patient_cell_grid %>%
  left_join(
    celltype_counts_all,
    by = c("patient_id", "response", "cell_type_analysis")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(patient_id, response) %>%
  mutate(
    total_cells = sum(n_cells),
    freq = n_cells / total_cells
  ) %>%
  ungroup()

write_csv(
  prop_all,
  "results/GSE120575_pre_celltype_proportions_all_immune_by_patient.csv"
)

prop_all

# =============================================================================
# Wilcoxon tests with BH correction: all immune cell proportions
# =============================================================================

stat_all <- prop_all %>%
  group_by(cell_type_analysis) %>%
  summarise(
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(freq[response == "Responder"], na.rm = TRUE),
    median_NR = median(freq[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = tryCatch(
      wilcox.test(freq ~ response, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

write_csv(
  stat_all,
  "results/GSE120575_pre_celltype_proportion_stats_all_immune_BH.csv"
)

stat_all

# =============================================================================
# Patient-level proportions within the T cell compartment
# =============================================================================

t_cell_types <- c(
  "Conventional / memory-like CD4 T",
  "Treg-like CD4 T",
  "CD4 T / Treg-enriched",
  "GZMK+ activated/transitional CD8 T",
  "TCF7+ memory-like CD8 T",
  "Cycling exhausted/tumor-reactive CD8 T",
  "gamma-delta / NK-like cytotoxic T"
)

t_cell_types_present <- t_cell_types[
  t_cell_types %in% unique(mel$cell_type_analysis)
]

meta_t <- mel@meta.data %>%
  filter(cell_type_analysis %in% t_cell_types_present)

patient_t_info <- meta_t %>%
  distinct(patient_id, response)

patient_t_grid <- patient_t_info %>%
  tidyr::crossing(cell_type_analysis = t_cell_types_present)

t_counts <- meta_t %>%
  count(patient_id, response, cell_type_analysis, name = "n_cells")

prop_t <- patient_t_grid %>%
  left_join(
    t_counts,
    by = c("patient_id", "response", "cell_type_analysis")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(patient_id, response) %>%
  mutate(
    total_t_cells = sum(n_cells),
    freq = n_cells / total_t_cells
  ) %>%
  ungroup()

write_csv(
  prop_t,
  "results/GSE120575_pre_Tcell_compartment_proportions_by_patient.csv"
)

prop_t

# =============================================================================
# Wilcoxon tests with BH correction: T cell compartment
# =============================================================================

stat_t <- prop_t %>%
  group_by(cell_type_analysis) %>%
  summarise(
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(freq[response == "Responder"], na.rm = TRUE),
    median_NR = median(freq[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = tryCatch(
      wilcox.test(freq ~ response, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

write_csv(
  stat_t,
  "results/GSE120575_pre_Tcell_compartment_proportion_stats_BH.csv"
)

stat_t

# =============================================================================
# Proportion plots
# =============================================================================

priority_cell_types <- c(
  "TCF7+ memory-like CD8 T",
  "Cycling exhausted/tumor-reactive CD8 T",
  "GZMK+ activated/transitional CD8 T",
  "Treg-like CD4 T",
  "Conventional / memory-like CD4 T",
  "CD4 T / Treg-enriched"
)

priority_cell_types_present <- priority_cell_types[
  priority_cell_types %in% unique(prop_all$cell_type_analysis)
]

p_prop_all_priority <- prop_all %>%
  filter(cell_type_analysis %in% priority_cell_types_present) %>%
  ggplot(aes(x = response, y = freq, fill = response)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_analysis, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Pre-treatment cell type proportions by response",
    subtitle = "Responder n = 9 patients · Non-responder n = 10 patients",
    x = NULL,
    y = "Proportion among all immune cells"
  )

p_prop_all_priority

ggsave(
  filename = "figures/GSE120575_pre_priority_celltype_proportions_all_immune.png",
  plot = p_prop_all_priority,
  width = 12,
  height = 7,
  dpi = 300
)

p_prop_t_priority <- prop_t %>%
  filter(cell_type_analysis %in% priority_cell_types_present) %>%
  ggplot(aes(x = response, y = freq, fill = response)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_analysis, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Pre-treatment T cell subset proportions by response",
    subtitle = "Responder n = 9 patients · Non-responder n = 10 patients",
    x = NULL,
    y = "Proportion within T cell compartment"
  )

p_prop_t_priority

ggsave(
  filename = "figures/GSE120575_pre_priority_Tcell_compartment_proportions.png",
  plot = p_prop_t_priority,
  width = 12,
  height = 7,
  dpi = 300
)

# 1) Cluster 0 subcluster sayısal kontrol
print(sub0_marker_positive_wide, n = Inf)

# 2) Tüm immune hücreler içinde oran istatistiği
print(stat_all, n = Inf)

# 3) Sadece T cell compartment içinde oran istatistiği
print(stat_t, n = Inf)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_analysis_ready.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_analysis_ready.rds")

# =============================================================================
# Denominator checks for proportion analyses
# =============================================================================

prop_all_check <- prop_all %>%
  group_by(patient_id, response) %>%
  summarise(
    sum_freq = sum(freq),
    n_cell_types = n(),
    .groups = "drop"
  )

prop_t_check <- prop_t %>%
  group_by(patient_id, response) %>%
  summarise(
    sum_freq = sum(freq),
    n_t_cell_types = n(),
    .groups = "drop"
  )

prop_all_check
prop_t_check

range(prop_all_check$sum_freq)
range(prop_t_check$sum_freq)

# =============================================================================
# Effect size summaries
# =============================================================================

stat_all_effect <- stat_all %>%
  mutate(
    ratio_R_over_NR = ifelse(median_NR > 0, median_R / median_NR, NA_real_),
    ratio_NR_over_R = ifelse(median_R > 0, median_NR / median_R, NA_real_)
  ) %>%
  arrange(p_value)

stat_t_effect <- stat_t %>%
  mutate(
    ratio_R_over_NR = ifelse(median_NR > 0, median_R / median_NR, NA_real_),
    ratio_NR_over_R = ifelse(median_R > 0, median_NR / median_R, NA_real_)
  ) %>%
  arrange(p_value)

print(stat_all_effect, n = Inf)
print(stat_t_effect, n = Inf)

write_csv(
  stat_all_effect,
  "results/GSE120575_pre_celltype_proportion_stats_all_immune_effects.csv"
)

write_csv(
  stat_t_effect,
  "results/GSE120575_pre_Tcell_compartment_proportion_stats_effects.csv"
)

# =============================================================================
# APC HLA-II module score
# =============================================================================

hla2_genes <- c(
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DPA1",
  "HLA-DPB1",
  "CD74"
)

hla2_genes_present <- hla2_genes[
  hla2_genes %in% rownames(mel)
]

hla2_genes_present

apc_cell_types <- c(
  "Myeloid / macrophage",
  "B cell",
  "pDC / dendritic-like"
)

apc <- subset(
  mel,
  subset = cell_type_analysis %in% apc_cell_types
)

apc <- AddModuleScore(
  apc,
  features = list(hla2_genes_present),
  name = "HLA2"
)

# AddModuleScore with name = "HLA2" creates column HLA21
head(apc@meta.data$HLA21)

# =============================================================================
# Patient-level APC HLA-II score
# =============================================================================

hla_apc_pt <- apc@meta.data %>%
  group_by(patient_id, response) %>%
  summarise(
    mean_hla2 = mean(HLA21, na.rm = TRUE),
    n_apc_cells = n(),
    .groups = "drop"
  )

hla_apc_pt

hla_apc_stat <- hla_apc_pt %>%
  summarise(
    comparison = "All APC",
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(mean_hla2[response == "Responder"], na.rm = TRUE),
    median_NR = median(mean_hla2[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = wilcox.test(mean_hla2 ~ response, exact = FALSE)$p.value
  )

hla_apc_stat

# =============================================================================
# Patient-level HLA-II score within myeloid / macrophage cells only
# =============================================================================

hla_myeloid_pt <- apc@meta.data %>%
  filter(cell_type_analysis == "Myeloid / macrophage") %>%
  group_by(patient_id, response) %>%
  summarise(
    mean_hla2 = mean(HLA21, na.rm = TRUE),
    n_myeloid_cells = n(),
    .groups = "drop"
  )

hla_myeloid_pt

hla_myeloid_stat <- hla_myeloid_pt %>%
  summarise(
    comparison = "Myeloid / macrophage",
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(mean_hla2[response == "Responder"], na.rm = TRUE),
    median_NR = median(mean_hla2[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = wilcox.test(mean_hla2 ~ response, exact = FALSE)$p.value
  )

hla_myeloid_stat

# =============================================================================
# Combine HLA-II tests and adjust p-values
# =============================================================================

hla_stats <- bind_rows(
  hla_apc_stat,
  hla_myeloid_stat
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  )

hla_stats

write_csv(
  hla_apc_pt,
  "results/GSE120575_pre_APC_HLAII_score_by_patient.csv"
)

write_csv(
  hla_myeloid_pt,
  "results/GSE120575_pre_myeloid_HLAII_score_by_patient.csv"
)

write_csv(
  hla_stats,
  "results/GSE120575_pre_HLAII_score_stats_BH.csv"
)

# =============================================================================
# HLA-II score plots
# =============================================================================

hla_plot_df <- bind_rows(
  hla_apc_pt %>%
    mutate(comparison = "All APC") %>%
    select(patient_id, response, comparison, mean_hla2),
  
  hla_myeloid_pt %>%
    mutate(comparison = "Myeloid / macrophage") %>%
    select(patient_id, response, comparison, mean_hla2)
)

p_hla2 <- ggplot(
  hla_plot_df,
  aes(x = response, y = mean_hla2, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ comparison, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Patient-level HLA-II module score in APC populations",
    subtitle = "Responder n = 9 patients · Non-responder n = 10 patients",
    x = NULL,
    y = "Mean HLA-II module score"
  )

p_hla2

ggsave(
  filename = "figures/GSE120575_pre_APC_HLAII_module_score_by_response.png",
  plot = p_hla2,
  width = 9,
  height = 5,
  dpi = 300
)

print(range(prop_all_check$sum_freq))
print(range(prop_t_check$sum_freq))

print(hla_stats, n = Inf)

file.exists("figures/GSE120575_pre_APC_HLAII_module_score_by_response.png")
file.exists("results/GSE120575_pre_HLAII_score_stats_BH.csv")
p_hla2

# =============================================================================
# APC composition and cell-type-specific HLA-II checks
# =============================================================================

apc_types <- c(
  "B cell",
  "Myeloid / macrophage",
  "pDC / dendritic-like"
)

apc_meta <- apc@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  filter(cell_type_analysis %in% apc_types)

# -----------------------------------------------------------------------------
# 1) APC composition per patient
# -----------------------------------------------------------------------------

apc_patient_info <- apc_meta %>%
  distinct(patient_id, response)

apc_counts <- apc_meta %>%
  count(patient_id, response, cell_type_analysis, name = "n_cells")

apc_composition <- apc_patient_info %>%
  tidyr::crossing(cell_type_analysis = apc_types) %>%
  left_join(
    apc_counts,
    by = c("patient_id", "response", "cell_type_analysis")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(patient_id, response) %>%
  mutate(
    total_apc_cells = sum(n_cells),
    apc_frac = n_cells / total_apc_cells
  ) %>%
  ungroup()

apc_composition_stats <- apc_composition %>%
  group_by(cell_type_analysis) %>%
  summarise(
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(apc_frac[response == "Responder"], na.rm = TRUE),
    median_NR = median(apc_frac[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = tryCatch(
      wilcox.test(apc_frac ~ response, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(apc_composition_stats, n = Inf)

write_csv(
  apc_composition,
  "results/GSE120575_pre_APC_composition_by_patient.csv"
)

write_csv(
  apc_composition_stats,
  "results/GSE120575_pre_APC_composition_stats_BH.csv"
)

# =============================================================================
# APC composition and cell-type-specific HLA-II checks - FIXED VERSION
# =============================================================================

apc_types <- c(
  "B cell",
  "Myeloid / macrophage",
  "pDC / dendritic-like"
)

# metadata içinde cell_id zaten olduğu için rownames'i cell_id diye eklemiyoruz.
# Bunun yerine barcode adını kullanıyoruz.
apc_meta <- tibble(
  barcode = rownames(apc@meta.data),
  patient_id = apc@meta.data[["patient_id"]],
  response = apc@meta.data[["response"]],
  cell_type_analysis = apc@meta.data[["cell_type_analysis"]],
  HLA21 = apc@meta.data[["HLA21"]]
) %>%
  filter(cell_type_analysis %in% apc_types)

# -----------------------------------------------------------------------------
# 1) APC composition per patient
# -----------------------------------------------------------------------------

apc_patient_info <- apc_meta %>%
  distinct(patient_id, response)

apc_counts <- apc_meta %>%
  count(patient_id, response, cell_type_analysis, name = "n_cells")

apc_composition <- apc_patient_info %>%
  tidyr::crossing(cell_type_analysis = apc_types) %>%
  left_join(
    apc_counts,
    by = c("patient_id", "response", "cell_type_analysis")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  group_by(patient_id, response) %>%
  mutate(
    total_apc_cells = sum(n_cells),
    apc_frac = n_cells / total_apc_cells
  ) %>%
  ungroup()

apc_composition_stats <- apc_composition %>%
  group_by(cell_type_analysis) %>%
  summarise(
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(apc_frac[response == "Responder"], na.rm = TRUE),
    median_NR = median(apc_frac[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = tryCatch(
      wilcox.test(apc_frac ~ response, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(apc_composition_stats, n = Inf)

write_csv(
  apc_composition,
  "results/GSE120575_pre_APC_composition_by_patient.csv"
)

write_csv(
  apc_composition_stats,
  "results/GSE120575_pre_APC_composition_stats_BH.csv"
)

# -----------------------------------------------------------------------------
# 2) HLA-II score within each APC type separately
# -----------------------------------------------------------------------------

hla_by_apc_type_pt <- apc_meta %>%
  group_by(patient_id, response, cell_type_analysis) %>%
  summarise(
    n_cells = n(),
    mean_hla2 = mean(HLA21, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5)

hla_by_apc_type_stats <- hla_by_apc_type_pt %>%
  group_by(cell_type_analysis) %>%
  summarise(
    n_R = n_distinct(patient_id[response == "Responder"]),
    n_NR = n_distinct(patient_id[response == "Non-responder"]),
    median_R = median(mean_hla2[response == "Responder"], na.rm = TRUE),
    median_NR = median(mean_hla2[response == "Non-responder"], na.rm = TRUE),
    diff_median_R_minus_NR = median_R - median_NR,
    p_value = if (n_distinct(response) == 2) {
      tryCatch(
        wilcox.test(mean_hla2 ~ response, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(hla_by_apc_type_stats, n = Inf)

write_csv(
  hla_by_apc_type_pt,
  "results/GSE120575_pre_HLAII_by_APC_type_patient_level.csv"
)

write_csv(
  hla_by_apc_type_stats,
  "results/GSE120575_pre_HLAII_by_APC_type_stats_BH.csv"
)

# =============================================================================
# APC composition plot
# =============================================================================

p_apc_composition <- ggplot(
  apc_composition,
  aes(x = response, y = apc_frac, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_analysis, scales = "free_y") +
  theme_bw() +
  labs(
    title = "APC compartment composition by response",
    subtitle = "Only patients contributing APCs are included",
    x = NULL,
    y = "Fraction within APC compartment"
  )

p_apc_composition

ggsave(
  filename = "figures/GSE120575_pre_APC_composition_by_response.png",
  plot = p_apc_composition,
  width = 10,
  height = 5,
  dpi = 300
)

# =============================================================================
# HLA-II score by APC type plot
# =============================================================================

p_hla_by_apc_type <- ggplot(
  hla_by_apc_type_pt,
  aes(x = response, y = mean_hla2, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_analysis, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Patient-level HLA-II module score by APC type",
    subtitle = "Patients with at least 5 cells of the given APC type are included",
    x = NULL,
    y = "Mean HLA-II module score"
  )

p_hla_by_apc_type

ggsave(
  filename = "figures/GSE120575_pre_HLAII_by_APC_type_response.png",
  plot = p_hla_by_apc_type,
  width = 10,
  height = 5,
  dpi = 300
)

print(apc_composition_stats, n = Inf)
print(hla_by_apc_type_stats, n = Inf)

file.exists("figures/GSE120575_pre_APC_composition_by_response.png")
file.exists("figures/GSE120575_pre_HLAII_by_APC_type_response.png")

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_final_analysis_ready.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_final_analysis_ready.rds")

# =============================================================================
# B cell HLA-II patient-level values
# =============================================================================

b <- subset(
  apc,
  subset = cell_type_analysis == "B cell"
)

b_pt <- b@meta.data %>%
  group_by(patient_id, response) %>%
  summarise(
    n_cells = n(),
    mean_hla2 = mean(HLA21, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  arrange(response, mean_hla2)

print(b_pt, n = Inf)

write_csv(
  b_pt,
  "results/GSE120575_pre_Bcell_HLAII_patient_values.csv"
)

p_bcell_hla2_patient <- ggplot(
  b_pt,
  aes(x = response, y = mean_hla2, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.9) +
  theme_bw() +
  labs(
    title = "B cell HLA-II module score by response",
    subtitle = "Patients with at least 5 B cells included",
    x = NULL,
    y = "Mean HLA-II module score"
  )

p_bcell_hla2_patient

ggsave(
  filename = "figures/GSE120575_pre_Bcell_HLAII_patient_values.png",
  plot = p_bcell_hla2_patient,
  width = 6,
  height = 5,
  dpi = 300
)

# =============================================================================
# Raw HLA-II genes in B cells: visual check
# =============================================================================

b_hla_genes <- c(
  "HLA-DRA",
  "CD74",
  "HLA-DPB1",
  "HLA-DPA1",
  "HLA-DRB1"
)

b_hla_genes_present <- b_hla_genes[
  b_hla_genes %in% rownames(b)
]

b_hla_genes_present

p_bcell_hla_genes_vln <- VlnPlot(
  b,
  features = b_hla_genes_present,
  group.by = "response",
  pt.size = 0,
  ncol = 3
)

p_bcell_hla_genes_vln

ggsave(
  filename = "figures/GSE120575_pre_Bcell_raw_HLAII_genes_by_response.png",
  plot = p_bcell_hla_genes_vln,
  width = 10,
  height = 6,
  dpi = 300
)

print(b_pt, n = Inf)

file.exists("results/GSE120575_pre_Bcell_HLAII_patient_values.csv")
file.exists("figures/GSE120575_pre_Bcell_HLAII_patient_values.png")
file.exists("figures/GSE120575_pre_Bcell_raw_HLAII_genes_by_response.png")

# =============================================================================
# Control 1: Is B cell HLA-II score associated with B cell count?
# =============================================================================

b <- subset(
  apc,
  subset = cell_type_analysis == "B cell"
)

b_pt <- b@meta.data %>%
  group_by(patient_id, response) %>%
  summarise(
    n_cells = n(),
    mean_hla2 = mean(HLA21, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  arrange(response, mean_hla2)

print(b_pt, n = Inf)

# Overall Spearman correlation
b_count_cor_overall <- cor.test(
  b_pt$n_cells,
  b_pt$mean_hla2,
  method = "spearman",
  exact = FALSE
)

b_count_cor_overall

# Within response groups
b_count_cor_by_response <- b_pt %>%
  group_by(response) %>%
  summarise(
    n_patients = n(),
    rho = cor(n_cells, mean_hla2, method = "spearman"),
    .groups = "drop"
  )

b_count_cor_by_response

# =============================================================================
# Control 2: Patient-level raw HLA-II gene expression in B cells
# =============================================================================

b_raw_genes <- c(
  "HLA-DRA",
  "CD74",
  "HLA-DPB1",
  "HLA-DPA1",
  "HLA-DRB1"
)

b_raw_genes_present <- b_raw_genes[
  b_raw_genes %in% rownames(b)
]

b_raw_genes_present

b_expr <- GetAssayData(
  b,
  assay = "RNA",
  layer = "data"
)[b_raw_genes_present, , drop = FALSE]

b_raw_cell <- as.data.frame(t(as.matrix(b_expr))) %>%
  rownames_to_column("cell_id") %>%
  left_join(
    b@meta.data %>%
      rownames_to_column("cell_id") %>%
      select(cell_id, patient_id, response),
    by = "cell_id"
  )

b_raw_pt <- b_raw_cell %>%
  group_by(patient_id, response) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(b_raw_genes_present),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  arrange(response)

print(b_raw_pt, n = Inf)

write_csv(
  b_raw_pt,
  "results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv"
)

# =============================================================================
# Control 2: Patient-level raw HLA-II gene expression in B cells - FIXED
# =============================================================================

b_raw_genes <- c(
  "HLA-DRA",
  "CD74",
  "HLA-DPB1",
  "HLA-DPA1",
  "HLA-DRB1"
)

b_raw_genes_present <- b_raw_genes[
  b_raw_genes %in% rownames(b)
]

b_raw_genes_present

b_expr <- GetAssayData(
  b,
  assay = "RNA",
  layer = "data"
)[b_raw_genes_present, , drop = FALSE]

# Expression matrix: cells x genes
b_expr_df <- as.data.frame(
  t(as.matrix(b_expr))
)

# Burada cell_id değil, barcode diyoruz ki metadata ile çakışmasın
b_expr_df$barcode <- rownames(b_expr_df)

b_meta_df <- tibble(
  barcode = rownames(b@meta.data),
  patient_id = b@meta.data[["patient_id"]],
  response = b@meta.data[["response"]]
)

b_raw_cell <- b_expr_df %>%
  left_join(
    b_meta_df,
    by = "barcode"
  )

b_raw_pt <- b_raw_cell %>%
  group_by(patient_id, response) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(b_raw_genes_present),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  arrange(response)

print(b_raw_pt, n = Inf)

write_csv(
  b_raw_pt,
  "results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv"
)

file.exists("results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv")

# =============================================================================
# Control 3: Which patients were included/excluded from B cell HLA-II analysis?
# =============================================================================

# All patients in the pre-treatment Seurat object
all_patient_info <- mel@meta.data %>%
  as_tibble() %>%
  distinct(patient_id, response)

# Number of B cells per patient
b_counts_all_patients <- mel@meta.data %>%
  as_tibble() %>%
  filter(cell_type_analysis == "B cell") %>%
  count(patient_id, response, name = "n_b_cells")

# Inclusion table for B cell HLA-II analysis
b_patient_inclusion <- all_patient_info %>%
  left_join(
    b_counts_all_patients,
    by = c("patient_id", "response")
  ) %>%
  mutate(
    n_b_cells = coalesce(n_b_cells, 0L),
    included_Bcell_HLAII = n_b_cells >= 5
  ) %>%
  arrange(response, included_Bcell_HLAII, n_b_cells)

print(b_patient_inclusion, n = Inf)

# Summary: how many patients were included/excluded per response group?
b_patient_inclusion_summary <- b_patient_inclusion %>%
  count(response, included_Bcell_HLAII, name = "n_patients") %>%
  arrange(response, included_Bcell_HLAII)

b_patient_inclusion_summary

# More explicit summary table
b_patient_inclusion_summary_wide <- b_patient_inclusion %>%
  group_by(response) %>%
  summarise(
    total_patients = n(),
    included_patients = sum(included_Bcell_HLAII),
    excluded_patients = sum(!included_Bcell_HLAII),
    median_b_cells = median(n_b_cells),
    min_b_cells = min(n_b_cells),
    max_b_cells = max(n_b_cells),
    .groups = "drop"
  )

b_patient_inclusion_summary_wide

# Save outputs
write_csv(
  b_patient_inclusion,
  "results/GSE120575_pre_Bcell_HLAII_patient_inclusion.csv"
)

write_csv(
  b_patient_inclusion_summary,
  "results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary.csv"
)

write_csv(
  b_patient_inclusion_summary_wide,
  "results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary_wide.csv"
)

# Plot B cell counts per patient by response
p_bcell_counts_inclusion <- ggplot(
  b_patient_inclusion,
  aes(x = response, y = n_b_cells, fill = included_Bcell_HLAII)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.9) +
  theme_bw() +
  labs(
    title = "B cell counts per patient and inclusion in B cell HLA-II analysis",
    subtitle = "Included if patient contributed at least 5 B cells",
    x = NULL,
    y = "Number of B cells",
    fill = "Included"
  )

p_bcell_counts_inclusion

ggsave(
  filename = "figures/GSE120575_pre_Bcell_HLAII_inclusion_counts.png",
  plot = p_bcell_counts_inclusion,
  width = 7,
  height = 5,
  dpi = 300
)

# File checks
file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion.csv")
file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary.csv")
file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary_wide.csv")
file.exists("figures/GSE120575_pre_Bcell_HLAII_inclusion_counts.png")

print(b_patient_inclusion, n = Inf)

b_patient_inclusion_summary_wide

file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion.csv")
file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary.csv")
file.exists("results/GSE120575_pre_Bcell_HLAII_patient_inclusion_summary_wide.csv")
file.exists("figures/GSE120575_pre_Bcell_HLAII_inclusion_counts.png")

# =============================================================================
# Save session information for reproducibility
# =============================================================================

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

file.exists("sessionInfo.txt")

# B cell count correlation
cor.test(b_pt$n_cells, b_pt$mean_hla2, method = "spearman", exact = FALSE)

b_pt %>%
  group_by(response) %>%
  summarise(
    n_patients = n(),
    rho = cor(n_cells, mean_hla2, method = "spearman"),
    .groups = "drop"
  )

# Check raw HLA-II patient-level file
file.exists("results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv")

# Save session info
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

file.exists("sessionInfo.txt")

# =============================================================================
# Final control: patient-level raw HLA-II gene expression in B cells
# =============================================================================

b_raw_genes <- c(
  "HLA-DRA",
  "CD74",
  "HLA-DPB1",
  "HLA-DRB1"
)

b_raw_genes_present <- b_raw_genes[
  b_raw_genes %in% rownames(b)
]

b_raw_genes_present

# Expression matrix: genes x cells
b_expr <- GetAssayData(
  b,
  assay = "RNA",
  layer = "data"
)[b_raw_genes_present, , drop = FALSE]

# Convert to cells x genes
b_expr_df <- as.data.frame(t(as.matrix(b_expr)))

# Avoid duplicate cell_id problem: use barcode
b_expr_df$barcode <- rownames(b_expr_df)

b_meta_df <- tibble(
  barcode = rownames(b@meta.data),
  patient_id = b@meta.data[["patient_id"]],
  response = b@meta.data[["response"]]
)

b_raw_cell <- b_expr_df %>%
  left_join(
    b_meta_df,
    by = "barcode"
  )

# Patient-level mean raw HLA-II gene expression
b_genes <- b_raw_cell %>%
  group_by(patient_id, response) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(b_raw_genes_present),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  arrange(response, HLA_DRA = `HLA-DRA`)

print(b_genes, n = Inf)

write_csv(
  b_genes,
  "results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv"
)

# =============================================================================
# Patient-level Wilcoxon tests for raw HLA-II genes in B cells
# =============================================================================

b_gene_stats <- map_dfr(
  b_raw_genes_present,
  function(g) {
    tibble(
      gene = g,
      n_R = n_distinct(b_genes$patient_id[b_genes$response == "Responder"]),
      n_NR = n_distinct(b_genes$patient_id[b_genes$response == "Non-responder"]),
      median_R = median(b_genes[[g]][b_genes$response == "Responder"], na.rm = TRUE),
      median_NR = median(b_genes[[g]][b_genes$response == "Non-responder"], na.rm = TRUE),
      diff_median_R_minus_NR = median_R - median_NR,
      p_value = wilcox.test(
        b_genes[[g]] ~ b_genes$response,
        exact = FALSE
      )$p.value
    )
  }
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(b_gene_stats, n = Inf)

write_csv(
  b_gene_stats,
  "results/GSE120575_pre_Bcell_raw_HLAII_gene_stats_BH.csv"
)

file.exists("results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv")
file.exists("results/GSE120575_pre_Bcell_raw_HLAII_gene_stats_BH.csv")

# =============================================================================
# Patient-level raw HLA-II gene expression plot in B cells
# =============================================================================

b_genes_long <- b_genes %>%
  pivot_longer(
    cols = all_of(b_raw_genes_present),
    names_to = "gene",
    values_to = "mean_expression"
  )

p_bcell_raw_hla_patient <- ggplot(
  b_genes_long,
  aes(x = response, y = mean_expression, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.3, alpha = 0.9) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Patient-level raw HLA-II gene expression in B cells",
    subtitle = "Patients with at least 5 B cells included",
    x = NULL,
    y = "Mean log-normalized expression"
  )

p_bcell_raw_hla_patient

ggsave(
  filename = "figures/GSE120575_pre_Bcell_raw_HLAII_patient_means.png",
  plot = p_bcell_raw_hla_patient,
  width = 10,
  height = 6,
  dpi = 300
)

file.exists("figures/GSE120575_pre_Bcell_raw_HLAII_patient_means.png")

print(b_genes, n = Inf)
print(b_gene_stats, n = Inf)

file.exists("results/GSE120575_pre_Bcell_raw_HLAII_patient_means.csv")
file.exists("results/GSE120575_pre_Bcell_raw_HLAII_gene_stats_BH.csv")
file.exists("figures/GSE120575_pre_Bcell_raw_HLAII_patient_means.png")

getwd()
list.files()
file.exists("README.md")
file.exists("sessionInfo.txt")
file.exists("figures")
file.exists("results")

list.files("figures")

main_figures <- c(
  "figures/01_umap_annotated.png" =
    "figures/GSE120575_pre_UMAP_celltype_final_annotation.png",
  
  "figures/02_dotplot_celltype_markers.png" =
    "figures/GSE120575_pre_annotation_marker_dotplot.png",
  
  "figures/03_tcell_proportions.png" =
    "figures/GSE120575_pre_priority_Tcell_compartment_proportions.png",
  
  "figures/04_apc_composition.png" =
    "figures/GSE120575_pre_APC_composition_by_response.png",
  
  "figures/05_apc_hla2_module_score.png" =
    "figures/GSE120575_pre_APC_HLAII_module_score_by_response.png",
  
  "figures/06_bcell_hla2_patient_level.png" =
    "figures/GSE120575_pre_Bcell_HLAII_patient_values.png"
)

file.copy(
  from = unname(main_figures),
  to = names(main_figures),
  overwrite = TRUE
)

file.exists(names(main_figures))

list.files("figures", pattern = "^[0-9]{2}_")

