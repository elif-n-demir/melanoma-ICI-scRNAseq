library(Seurat)
library(tidyverse)
library(Matrix)

metadata_file <- "results/GSE120575_melanoma_metadata_clean.csv"

expr_file_txt <- "data_raw/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"

melanoma_metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

file.exists(expr_file_txt)

con <- file(expr_file_txt, open = "r")

header_line <- readLines(con, n = 1)
timepoint_line <- readLines(con, n = 1)

close(con)

expr_cell_ids_all <- strsplit(
  header_line,
  "\t",
  fixed = TRUE
)[[1]][-1]

expr_timepoints_all <- strsplit(
  timepoint_line,
  "\t",
  fixed = TRUE
)[[1]][-1]

length(expr_cell_ids_all)

head(expr_cell_ids_all)

sum(expr_cell_ids_all %in% melanoma_metadata$cell_id)

all(expr_cell_ids_all %in% melanoma_metadata$cell_id)

melanoma_metadata %>%
  count(timepoint, response)

melanoma_metadata %>%
  count(therapy, response)

melanoma_metadata %>%
  count(patient_id, timepoint, response) %>%
  arrange(patient_id, timepoint)

melanoma_metadata %>%
  count(timepoint, response) %>%
  print(n = Inf)

melanoma_metadata %>%
  count(therapy, timepoint, response) %>%
  print(n = Inf)

set.seed(123)

selected_metadata <- melanoma_metadata %>%
  filter(
    timepoint == "Pre",
    response %in% c("Responder", "Non-responder")
  ) %>%
  group_by(response) %>%
  slice_sample(n = 1000) %>%
  ungroup()

selected_cells <- selected_metadata$cell_id

selected_col_positions <- match(
  selected_cells,
  expr_cell_ids_all
) + 1

length(selected_cells)

selected_metadata %>%
  count(response)

sum(is.na(selected_col_positions))

head(selected_cells)

head(selected_col_positions)

selected_metadata <- melanoma_metadata %>%
  filter(timepoint == "Pre")

selected_cells <- selected_metadata$cell_id

selected_col_positions <- match(
  selected_cells,
  expr_cell_ids_all
) + 1

nrow(selected_metadata)

selected_metadata %>%
  count(response)

selected_metadata %>%
  count(patient_id, response) %>%
  arrange(patient_id)

sum(is.na(selected_col_positions))

head(selected_cells)

head(selected_col_positions)

# Alignment test before reading the full selected expression matrix

con <- file(expr_file_txt, "r")

# Skip first two lines:
# 1 = cell IDs
# 2 = timepoint row
invisible(readLines(con, n = 2))

# Read first gene row only
test_line <- readLines(con, n = 1)

close(con)

parts <- strsplit(
  test_line,
  "\t",
  fixed = TRUE
)[[1]]

length(parts)

parts[1]

parts[selected_col_positions[1:5]]

# Read selected Pre-treatment expression matrix safely

con <- file(expr_file_txt, open = "r")

# Skip first two lines:
# 1 = cell IDs
# 2 = timepoint row
invisible(readLines(con, n = 2))

sparse_chunks <- list()
gene_names <- character()

chunk_size <- 2000
total_read <- 0

repeat {
  lines <- readLines(con, n = chunk_size)
  
  if (length(lines) == 0) {
    break
  }
  
  parsed <- strsplit(
    lines,
    "\t",
    fixed = TRUE
  )
  
  genes <- vapply(
    parsed,
    function(x) x[1],
    character(1)
  )
  
  vals <- t(vapply(
    parsed,
    function(x) as.numeric(x[selected_col_positions]),
    numeric(length(selected_col_positions))
  ))
  
  # Remove genes with zero expression across all selected Pre-treatment cells
  keep <- rowSums(vals) > 0
  
  if (any(keep)) {
    vals_sparse <- Matrix(
      vals[keep, , drop = FALSE],
      sparse = TRUE
    )
    
    rownames(vals_sparse) <- genes[keep]
    
    sparse_chunks[[length(sparse_chunks) + 1]] <- vals_sparse
    gene_names <- c(gene_names, genes[keep])
  }
  
  total_read <- total_read + length(lines)
  
  cat(
    "Okunan gen:", total_read,
    "| tutulan gen:", length(gene_names),
    "\n"
  )
}

close(con)

expr_sparse <- do.call(
  rbind,
  sparse_chunks
)

rownames(expr_sparse) <- make.unique(rownames(expr_sparse))
colnames(expr_sparse) <- selected_cells

dim(expr_sparse)

expr_sparse[1:5, 1:5]

saveRDS(
  expr_sparse,
  "data_processed/GSE120575_pre_expr_sparse.rds"
)

write_csv(
  selected_metadata,
  "results/GSE120575_pre_metadata_selected.csv"
)

rm(sparse_chunks, vals, vals_sparse, parsed, lines)
gc()

dim(expr_sparse)

nrow(expr_sparse)
ncol(expr_sparse)

file.exists("data_processed/GSE120575_pre_expr_sparse.rds")
file.exists("results/GSE120575_pre_metadata_selected.csv")

meta_df <- as.data.frame(selected_metadata)

rownames(meta_df) <- meta_df$cell_id

meta_df <- meta_df[colnames(expr_sparse), ]

stopifnot(
  all(rownames(meta_df) == colnames(expr_sparse))
)

mel <- CreateSeuratObject(
  counts = expr_sparse,
  meta.data = meta_df,
  project = "melanoma_ICI_pre",
  min.cells = 3
)

mel

table(mel$response)

table(mel$patient_id)

saveRDS(
  mel,
  "data_processed/GSE120575_pre_seurat_raw.rds"
)

file.exists("data_processed/GSE120575_pre_seurat_raw.rds")