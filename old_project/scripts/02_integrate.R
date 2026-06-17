#!/usr/bin/env Rscript
# 02_integrate.R
# fastMNN integration (temporal merge order) + timepoint mixing-matrix diagnostic.
# Produces both baseline (PCA) and post-MNN matrices so integration quality
# is checked every run — not just once during exploration.
# Input : 01_merged.rds
# Output: 02_integrated.rds | 02_mixing_*.csv | 02_mixing_comparison.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(batchelor)
  library(BiocParallel)
  library(BiocNeighbors)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",  type = "character"),
  make_option("--n_dims", type = "integer", default = 30L)
)))

tp_order <- c("E80", "E825", "E95", "E105")
tp_cols  <- setNames(viridisLite::viridis(4, option = "D", direction = -1), tp_order)

obj <- readRDS(opt$input)
hvg <- VariableFeatures(obj)

# ── Helper: timepoint mixing matrix ──────────────────────────────────────────
# For each focal timepoint, average the neighbour-label distribution (k = 30).
# Diagonal = self-containment; off-diagonal = cross-timepoint mixing.
mixing_matrix <- function(emb, labels, k = 30) {
  nn  <- findKNN(emb, k = k)$index
  map_dfr(tp_order, function(tp) {
    cells <- which(labels == tp)
    nbrs  <- as.character(labels[nn[cells, ]])
    tibble(focal = tp, neighbour = tp_order,
           frac  = map_dbl(tp_order, ~ mean(nbrs == .x)))
  })
}

# ── Baseline: PCA (no integration) ───────────────────────────────────────────
obj <- ScaleData(obj, verbose = FALSE) |>
  RunPCA(npcs = opt$n_dims, verbose = FALSE)

mix_base <- mixing_matrix(Embeddings(obj, "pca")[, 1:opt$n_dims],
                          as.character(obj$timepoint))
write_csv(mix_base, "02_mixing_baseline.csv")

# ── fastMNN integration ───────────────────────────────────────────────────────
# Temporal merge order preserves the developmental gradient:
# E80 → E825 → E95 → E105 (do not reorder).
log_mat <- GetAssayData(obj, layer = "data")
mnn_out <- fastMNN(log_mat, batch = obj$timepoint,
                   subset.row = hvg, d = opt$n_dims,
                   merge.order = tp_order,
                   BPPARAM = SerialParam())   # serial avoids Mac ARM/Rosetta issues

corrected           <- reducedDim(mnn_out, "corrected")
rownames(corrected) <- colnames(obj)
obj[["mnn"]] <- CreateDimReducObject(corrected, key = "MNN_")
obj <- RunUMAP(obj, reduction = "mnn", dims = seq_len(opt$n_dims),
               reduction.name = "umap.mnn", seed.use = 42, verbose = FALSE)

mix_mnn <- mixing_matrix(corrected, as.character(obj$timepoint))
write_csv(mix_mnn, "02_mixing_mnn.csv")

# ── Comparison plot ───────────────────────────────────────────────────────────
p <- bind_rows(
  mutate(mix_base, arm = "Baseline (no integration)"),
  mutate(mix_mnn,  arm = "fastMNN")
) |>
  mutate(across(c(focal, neighbour), ~ factor(.x, tp_order))) |>
  ggplot(aes(neighbour, focal, fill = frac, label = round(frac, 2))) +
  geom_tile() + geom_text(size = 3.5, color = "white") +
  scale_fill_viridis_c(option = "D") +
  facet_wrap(~arm) +
  labs(title = "Timepoint mixing matrix — baseline vs fastMNN",
       subtitle = "Diagonal = self-containment | off-diagonal = cross-timepoint mixing",
       x = "neighbour timepoint", y = "focal timepoint") +
  theme_minimal(base_size = 11)
ggsave("02_mixing_comparison.pdf", p, width = 12, height = 5)

saveRDS(obj, "02_integrated.rds")
message("[integrate] done → 02_integrated.rds")
