#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 07: DIM_REDUCTION_AND_CLUSTER
# -----------------------------------------------------------------------------
# PCA -> UMAP -> Louvain clustering on the SCT-normalized object, plus an elbow
# plot (PC-count sanity) and a UMAP (cluster sanity).
#
# PLACE IN PIPELINE:
#   06 NORMALIZE_DATA -> [07 DIM_REDUCTION_AND_CLUSTER]   (last step of Phase 1)
#
# INPUT  : 06_seurat_normalized.rds
# OUTPUT : 07_seurat_clustered.rds   (-> Phase 2 annotation)
#          07_elbow.pdf  (PC variance, with the chosen n_pcs marked)
#          07_umap.pdf   (clusters)
# ENV    : paa-dev.yml  (Seurat 5.x, tidyverse, fs, glue, optparse)
#
# IMPORTANT — these clusters are a QC DIAGNOSTIC, not an analytical deliverable.
#   The clustering that feeds results happens on the integrated atlas (Phase 3)
#   and the caudal-arch subset (Phase 4a). A single sensible resolution is fine
#   here; this is why clustree was dropped from Phase 1.
#
# REPRODUCIBILITY: every stochastic call (PCA, UMAP, clustering) is seeded.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(tidyverse)   # ggplot2 for the diagnostic plots
  library(fs)
  library(glue)
})

# ---- 2. Command-line arguments ---------------------------------------------
option_list <- list(
  make_option("--in_rds",      type = "character", help = "06_seurat_normalized.rds"),
  make_option("--n_pcs",       type = "integer", default = 30L,
              help = "PCs used for UMAP + graph [default %default] (validate via elbow)"),
  make_option("--cluster_res", type = "double", default = 0.4,
              help = "Louvain resolution [default %default]"),
  make_option("--seed",        type = "integer", default = 42L,
              help = "RNG seed for PCA/UMAP/clustering [default %default]"),
  make_option("--out_rds",     type = "character", help = "Output 07_seurat_clustered.rds"),
  make_option("--out_elbow",   type = "character", help = "Output 07_elbow.pdf"),
  make_option("--out_umap",    type = "character", help = "Output 07_umap.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("in_rds", "out_rds", "out_elbow", "out_umap")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))
if (!file_exists(opt$in_rds)) stop(glue("Input not found: {opt$in_rds}"))

# ---- 3. Load ----------------------------------------------------------------
message(glue("[07] Loading {opt$in_rds}"))
obj <- readRDS(opt$in_rds)
sid <- obj$sample_id[1]

# ---- 4. PCA (compute 50; use n_pcs downstream) ------------------------------
set.seed(opt$seed)
obj <- RunPCA(obj, npcs = 50, seed.use = opt$seed, verbose = FALSE)

# Elbow: where added PCs stop explaining meaningful variance. The dashed line is
# the chosen n_pcs — inspect that it sits at/after the elbow.
p_elbow <- ElbowPlot(obj, ndims = 50) +
  geom_vline(xintercept = opt$n_pcs, linetype = "dashed", colour = "firebrick") +
  labs(title = glue("{sid} — PCA elbow (n_pcs = {opt$n_pcs})")) +
  theme_bw()
dir_create(path_dir(opt$out_elbow))
ggsave(opt$out_elbow, p_elbow, width = 6, height = 4.5)

# ---- 5. UMAP + Louvain clustering (all seeded) ------------------------------
dims <- 1:opt$n_pcs
obj <- RunUMAP(obj, dims = dims, seed.use = opt$seed, verbose = FALSE)
obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = opt$cluster_res, random.seed = opt$seed, verbose = FALSE)

p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(title = glue("{sid} — clusters (res {opt$cluster_res}, {opt$n_pcs} PCs)")) +
  theme_bw()
dir_create(path_dir(opt$out_umap))
ggsave(opt$out_umap, p_umap, width = 7, height = 6)

# ---- 6. Report --------------------------------------------------------------
n_clust <- length(levels(obj$seurat_clusters))
message(glue(
  "[07] '{sid}': {n_clust} clusters at resolution {opt$cluster_res} ",
  "(n_pcs = {opt$n_pcs}, {ncol(obj)} cells). Plots: ",
  "{path_file(opt$out_elbow)}, {path_file(opt$out_umap)}"
))

# ---- 7. Save ----------------------------------------------------------------
saveRDS(obj, opt$out_rds)
message(glue("[07] Wrote {opt$out_rds}"))
