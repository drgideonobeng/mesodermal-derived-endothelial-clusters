#!/usr/bin/env Rscript
# 03_cluster.R
# Build kNN graph on the MNN-corrected embedding and cluster at a fixed
# resolution. Resolution is a pipeline parameter — change it in pipeline.yml.
# To re-run the interactive resolution sweep (clustree), use the --run_sweep
# flag: the sweep runs but does NOT replace the locked resolution.
# Input : 02_integrated.rds
# Output: 03_clustered.rds | 03_umap_clusters.pdf | 03_cluster_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(clustree)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type = "character"),
  make_option("--n_dims",     type = "integer", default = 30L),
  make_option("--resolution", type = "double",  default = 0.4),
  make_option("--algorithm",  type = "integer", default = 1L),
  make_option("--run_sweep",  action = "store_true", default = FALSE,
              help = "Produce a clustree plot across resolutions 0.1-1.2")
)))

set.seed(42)
tp_order <- c("E80", "E825", "E95", "E105")
tp_cols  <- setNames(viridisLite::viridis(4, option = "D", direction = -1), tp_order)

obj  <- readRDS(opt$input)
dims <- seq_len(opt$n_dims)

# ── Graph construction ────────────────────────────────────────────────────────
obj <- FindNeighbors(obj, reduction = "mnn", dims = dims,
                     graph.name = c("mnn_nn", "mnn_snn"), verbose = FALSE)

# ── Optional: resolution sweep for exploratory use ───────────────────────────
if (opt$run_sweep) {
  resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2)
  for (res in resolutions) {
    obj <- FindClusters(obj, graph.name = "mnn_snn",
                        resolution = res, algorithm = opt$algorithm, verbose = FALSE)
    obj[[glue("clust_res{res}")]] <- Idents(obj)
  }
  ct <- clustree(obj@meta.data, prefix = "clust_res") + ggtitle("Resolution sweep")
  ggsave("03_clustree_sweep.pdf", ct, width = 10, height = 14)
  message("[cluster] clustree sweep written → 03_clustree_sweep.pdf")
}

# ── Cluster at locked resolution ──────────────────────────────────────────────
res_col <- glue("clust_res{opt$resolution}")
obj <- FindClusters(obj, graph.name = "mnn_snn",
                    resolution = opt$resolution, algorithm = opt$algorithm,
                    verbose = FALSE)
obj[[res_col]] <- Idents(obj)
message(glue("[cluster] res {opt$resolution} → {nlevels(Idents(obj))} clusters"))

# ── UMAP plots ────────────────────────────────────────────────────────────────
p_clust <- DimPlot(obj, reduction = "umap.mnn", label = TRUE, repel = TRUE) +
  ggtitle(glue("Clusters — res {opt$resolution}")) + NoLegend()
p_time  <- DimPlot(obj, reduction = "umap.mnn",
                   group.by = "timepoint", cols = tp_cols) +
  ggtitle("Timepoint")
ggsave("03_umap_clusters.pdf", p_clust + p_time, width = 14, height = 6)

# ── Cluster summary ───────────────────────────────────────────────────────────
obj@meta.data |>
  as_tibble() |>
  count(cluster = .data[[res_col]], timepoint) |>
  write_csv("03_cluster_summary.csv")

saveRDS(obj, "03_clustered.rds")
message("[cluster] done → 03_clustered.rds")
