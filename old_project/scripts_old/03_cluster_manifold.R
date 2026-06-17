#!/usr/bin/env Rscript
# ============================================================================
# 03_cluster_manifold.R
# ----------------------------------------------------------------------------
# STEP 1 of whole-manifold annotation: define robust clusters on the fastMNN
# embedding before we attempt to label anything.
#
# Strategy (why this, not a single FindClusters call):
#   A single resolution is an arbitrary choice. Instead we SWEEP a range of
#   resolutions, store each as its own metadata column, and use clustree to
#   see how clusters split / merge as resolution increases. We pick a
#   resolution on a STABLE PLATEAU -- where raising the resolution only
#   cleanly subdivides already-coherent clusters, rather than reshuffling
#   cells across parents (which shows up as crossing edges in the tree).
#   The SC3 stability index, overlaid on a second clustree panel, gives a
#   per-cluster robustness score across the sweep.
#
# Inputs : results/objects/02_merged_mnn.rds   (Seurat obj w/ 'mnn' reduction)
# Outputs: results/objects/03_clustered.rds     (all resolution cols stored)
#          results/plots/03_clustree.pdf         (sweep tree, x2 panels)
#          results/plots/03_umap_grid_resolutions.pdf
#          results/plots/03_umap_clusters_preview.pdf
#          results/plots/03_cluster_timepoint_composition.pdf
#          results/03_cluster_counts_by_resolution.csv
#          results/03_cluster_timepoint_composition.csv
#
# Usage (in the paa-dev env, from the project root):
#   conda run -n paa-dev Rscript scripts/03_cluster_manifold.R \
#       --input  results/objects/02_merged_mnn.rds \
#       --outdir results
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clustree)
  library(patchwork)
  library(glue)
  library(optparse)
})

# ---- CLI -------------------------------------------------------------------
option_list <- list(
  make_option("--input",     type = "character",
              default = "results/objects/02_merged_mnn.rds",
              help = "fastMNN-integrated Seurat object [%default]"),
  make_option("--outdir",    type = "character", default = "results",
              help = "Output root (expects objects/ and plots/ subdirs) [%default]"),
  make_option("--reduction", type = "character", default = "mnn",
              help = "Reduction to build the neighbour graph on [%default]"),
  make_option("--umap",      type = "character", default = "umap.mnn",
              help = "UMAP reduction for plotting [%default]"),
  make_option("--dims",      type = "integer",   default = 30L,
              help = "Number of mnn dimensions to use [%default]"),
  make_option("--algorithm", type = "integer",   default = 1L,
              help = "FindClusters algo: 1=Louvain, 2=Louvain+refine, 4=Leiden [%default]"),
  make_option("--resolutions", type = "character",
              default = "0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2",
              help = "Comma-separated resolution sweep [%default]"),
  make_option("--resolution_preview", type = "double", default = 0.4,
              help = "Provisional resolution rendered as the cluster UMAP [%default]"),
  make_option("--timepoint_col", type = "character", default = "timepoint",
              help = "Metadata column holding the developmental stage [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir_obj  <- file.path(opt$outdir, "objects")
dir_plot <- file.path(opt$outdir, "plots")
dir.create(dir_obj,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)

resolutions <- as.numeric(strsplit(opt$resolutions, ",")[[1]])
tp_order <- c("E80", "E825", "E95", "E105")
tp_cols  <- setNames(viridisLite::viridis(4, option = "D", direction = -1), tp_order)

set.seed(42)  # FindClusters has a stochastic component; pin it

# ---- Load ------------------------------------------------------------------
message(glue("[load] {opt$input}"))
obj <- readRDS(opt$input)

stopifnot(opt$reduction %in% Reductions(obj))
if (!opt$timepoint_col %in% colnames(obj@meta.data)) {
  stop(glue("Timepoint column '{opt$timepoint_col}' not found in meta.data. ",
            "Available: {paste(colnames(obj@meta.data), collapse=', ')}"))
}
# normalise stage levels for ordered plotting
obj@meta.data[[opt$timepoint_col]] <-
  factor(as.character(obj@meta.data[[opt$timepoint_col]]), levels = tp_order)

n_mnn <- ncol(Embeddings(obj, opt$reduction))
dims_use <- seq_len(min(opt$dims, n_mnn))
message(glue("[graph] FindNeighbors on '{opt$reduction}' dims 1:{max(dims_use)} ",
             "({n_mnn} available)"))

# ---- Neighbour graph on the corrected (mnn) space --------------------------
obj <- FindNeighbors(
  obj,
  reduction  = opt$reduction,
  dims       = dims_use,
  graph.name = c("mnn_nn", "mnn_snn"),
  verbose    = FALSE
)

# ---- Resolution sweep ------------------------------------------------------
# Each resolution is stored under a clean, clustree-friendly prefix.
res_cols <- character(0)
for (res in resolutions) {
  obj <- FindClusters(
    obj,
    graph.name = "mnn_snn",
    resolution = res,
    algorithm  = opt$algorithm,
    verbose    = FALSE
  )
  col <- glue("clust_res{res}")
  obj[[col]] <- Idents(obj)
  res_cols <- c(res_cols, col)
  message(glue("[sweep] res {res} -> {nlevels(Idents(obj))} clusters"))
}

# ---- n_clusters vs resolution ---------------------------------------------
counts_tbl <- tibble(
  resolution = resolutions,
  n_clusters = map_int(res_cols, ~ nlevels(obj@meta.data[[.x]]))
)
write_csv(counts_tbl, file.path(opt$outdir, "03_cluster_counts_by_resolution.csv"))
message("[table] cluster counts by resolution:")
print(counts_tbl)

# ---- clustree: structure + SC3 stability -----------------------------------
# Panel 1: node size = n cells, coloured by resolution -> shows split/merge.
# Panel 2: nodes coloured by SC3 stability -> high = robust across the sweep.
ct_struct <- clustree(obj@meta.data, prefix = "clust_res") +
  ggtitle("Cluster structure across resolutions")
ct_stab <- clustree(obj@meta.data, prefix = "clust_res",
                    node_colour = "sc3_stability") +
  ggtitle("SC3 stability (higher = more robust)")

ggsave(file.path(dir_plot, "03_clustree.pdf"),
       ct_struct / ct_stab, width = 11, height = 16)
message("[plot] 03_clustree.pdf")

# ---- UMAP grid across all resolutions --------------------------------------
grid_plots <- map2(res_cols, resolutions, function(col, res) {
  DimPlot(obj, reduction = opt$umap, group.by = col, label = TRUE,
          label.size = 3, repel = TRUE) +
    NoLegend() +
    ggtitle(glue("res {res}  |  k={nlevels(obj@meta.data[[col]])}")) +
    theme(plot.title = element_text(size = 9))
})
ggsave(file.path(dir_plot, "03_umap_grid_resolutions.pdf"),
       wrap_plots(grid_plots, ncol = 3),
       width = 15, height = 5 * ceiling(length(res_cols) / 3))
message("[plot] 03_umap_grid_resolutions.pdf")

# ---- Provisional resolution preview ----------------------------------------
preview_col <- glue("clust_res{opt$resolution_preview}")
if (!preview_col %in% res_cols) {
  warning(glue("Preview resolution {opt$resolution_preview} not in sweep; ",
               "falling back to {res_cols[1]}"))
  preview_col <- res_cols[1]
}
Idents(obj) <- preview_col
message(glue("[preview] provisional identity = {preview_col} ",
             "({nlevels(Idents(obj))} clusters) -- REVISABLE after clustree"))

p_clusters <- DimPlot(obj, reduction = opt$umap, label = TRUE, repel = TRUE) +
  ggtitle(glue("Clusters @ {preview_col} (provisional)"))
ggsave(file.path(dir_plot, "03_umap_clusters_preview.pdf"),
       p_clusters, width = 8, height = 7)
message("[plot] 03_umap_clusters_preview.pdf")

# ---- Cluster x timepoint composition ---------------------------------------
# In a timecourse this is diagnostic: stage-restricted clusters (one bar
# dominated by a single stage) are late- or early-emerging types; shared
# clusters spanning stages are the persistent backbone of the trajectory.
comp <- obj@meta.data |>
  as_tibble() |>
  count(cluster = .data[[preview_col]], stage = .data[[opt$timepoint_col]],
        name = "n") |>
  group_by(cluster) |>
  mutate(frac = n / sum(n),
         cluster_size = sum(n)) |>
  ungroup()
write_csv(comp, file.path(opt$outdir, "03_cluster_timepoint_composition.csv"))

p_comp <- ggplot(comp, aes(x = cluster, y = frac, fill = stage)) +
  geom_col() +
  scale_fill_manual(values = tp_cols, drop = FALSE) +
  labs(x = glue("cluster ({preview_col})"), y = "fraction of cluster",
       fill = "stage",
       title = "Cluster composition by developmental stage",
       subtitle = "Single-stage bars = stage-restricted types; mixed bars = persistent") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(dir_plot, "03_cluster_timepoint_composition.pdf"),
       p_comp, width = 10, height = 5)
message("[plot] 03_cluster_timepoint_composition.pdf")

# ---- Save ------------------------------------------------------------------
saveRDS(obj, file.path(dir_obj, "03_clustered.rds"))
message(glue("[save] {file.path(dir_obj, '03_clustered.rds')}"))
message("[done] Step 1 complete. Inspect 03_clustree.pdf to choose a resolution, ",
        "then we move to Step 2 (cluster DE + marker scoring).")