# Melanoma spatial transcriptomics demo
# Independent public human melanoma Visium dataset
# This is a workflow demonstration, not spatial validation of the GSE120575 responder/non-responder analysis.

library(Seurat)
library(ggplot2)

# -----------------------------
# 1. Paths
# -----------------------------

base_dir <- "extensions/melanoma_spatial_transcriptomics_demo"

data_dir <- file.path(
  base_dir,
  "data_raw",
  "10x_human_melanoma_if_ffpe"
)

fig_dir <- file.path(base_dir, "figures")
res_dir <- file.path(base_dir, "results")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

h5_file <- "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5"

# -----------------------------
# 2. Load spatial data
# -----------------------------

mel_spatial <- Load10X_Spatial(
  data.dir = data_dir,
  filename = h5_file
)

DefaultAssay(mel_spatial) <- "Spatial"

# Save basic object summary
object_summary <- data.frame(
  n_features = nrow(mel_spatial),
  n_spots = ncol(mel_spatial),
  assay = DefaultAssay(mel_spatial)
)

write.csv(
  object_summary,
  file = file.path(res_dir, "01_object_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# 3. Basic QC metadata
# -----------------------------

mel_spatial[["percent.mt"]] <- PercentageFeatureSet(
  mel_spatial,
  pattern = "^MT-"
)

qc_table <- mel_spatial@meta.data

write.csv(
  qc_table,
  file = file.path(res_dir, "02_spot_qc_metadata.csv")
)

# QC plots
qc_plot <- VlnPlot(
  mel_spatial,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
)

ggsave(
  filename = file.path(fig_dir, "01_qc_violin_plots.png"),
  plot = qc_plot,
  width = 11,
  height = 4,
  dpi = 300
)

# -----------------------------
# 4. Standard normalisation and clustering
# -----------------------------
# SCTransform was not used in this demo to avoid model-fitting warnings
# and to keep the workflow simple and reproducible.

mel_spatial <- NormalizeData(
  mel_spatial,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

mel_spatial <- FindVariableFeatures(
  mel_spatial,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

mel_spatial <- ScaleData(
  mel_spatial,
  verbose = FALSE
)

mel_spatial <- RunPCA(
  mel_spatial,
  npcs = 30,
  verbose = FALSE
)

mel_spatial <- FindNeighbors(
  mel_spatial,
  dims = 1:20,
  verbose = FALSE
)

mel_spatial <- FindClusters(
  mel_spatial,
  resolution = 0.4,
  verbose = FALSE
)

mel_spatial <- RunUMAP(
  mel_spatial,
  dims = 1:20,
  verbose = FALSE
)

# -----------------------------
# 5. Spatial cluster plot
# -----------------------------

cluster_plot <- SpatialDimPlot(
  mel_spatial,
  label = TRUE,
  label.size = 3
) +
  ggtitle("Spatial clustering of independent human melanoma Visium sample")

ggsave(
  filename = file.path(fig_dir, "02_spatial_clusters.png"),
  plot = cluster_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# -----------------------------
# 6. Marker spatial feature plots
# -----------------------------

marker_genes <- c(
  "SOX10", "MLANA", "PMEL", "S100B",
  "CD3D", "CD3E", "CD8A",
  "CD74", "HLA-DRA",
  "S100A8", "S100A9"
)

marker_genes_present <- intersect(marker_genes, rownames(mel_spatial))

write.csv(
  data.frame(marker = marker_genes, present = marker_genes %in% rownames(mel_spatial)),
  file = file.path(res_dir, "03_marker_gene_presence.csv"),
  row.names = FALSE
)

marker_plot <- SpatialFeaturePlot(
  mel_spatial,
  features = marker_genes_present,
  ncol = 3
)

ggsave(
  filename = file.path(fig_dir, "03_spatial_marker_featureplots.png"),
  plot = marker_plot,
  width = 13,
  height = 10,
  dpi = 300
)

# -----------------------------
# 7. Targeted antigen-presentation module score
# -----------------------------

ap_genes <- c(
  "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
  "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA",
  "HLA-DMA", "HLA-DMB"
)

ap_genes_present <- intersect(ap_genes, rownames(mel_spatial))

write.csv(
  data.frame(gene = ap_genes, present = ap_genes %in% rownames(mel_spatial)),
  file = file.path(res_dir, "04_antigen_presentation_gene_presence.csv"),
  row.names = FALSE
)

if (length(ap_genes_present) >= 2) {
  mel_spatial <- AddModuleScore(
    object = mel_spatial,
    features = list(ap_genes_present),
    assay = "Spatial",
    name = "AP_HLAII_Score"
  )
  
  ap_score_plot <- SpatialFeaturePlot(
    mel_spatial,
    features = "AP_HLAII_Score1"
  ) +
    ggtitle("Targeted antigen-presentation module score")
  
  ggsave(
    filename = file.path(fig_dir, "04_antigen_presentation_module_score.png"),
    plot = ap_score_plot,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  ap_score_table <- mel_spatial@meta.data[, c("AP_HLAII_Score1", "seurat_clusters")]
  
  write.csv(
    ap_score_table,
    file = file.path(res_dir, "05_antigen_presentation_module_score_by_spot.csv")
  )
}

# -----------------------------
# 8. Save processed object metadata only
# -----------------------------
# Raw data and full Seurat objects are not committed to GitHub.
# Only lightweight result tables and figures are saved.

metadata_out <- mel_spatial@meta.data

write.csv(
  metadata_out,
  file = file.path(res_dir, "06_processed_spot_metadata.csv")
)

# -----------------------------
# 9. Session information
# -----------------------------

sink(file.path(res_dir, "07_sessionInfo.txt"))
sessionInfo()
sink()

message("Spatial transcriptomics demo completed successfully.")
