#!/usr/bin/env Rscript
# 05_markers.R
# Unbiased cluster differential expression (Wilcoxon, positive markers only).
# Isolated from annotation so the slow DE step can be cached independently.
# Input : 04_cc_scored.rds
# Output: 05_markers_all.csv | 05_markers_top.csv | 05_dotplot_markers.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(glue)
  library(optparse)
  library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type = "character"),
  make_option("--resolution", type = "double",  default = 0.4),
  make_option("--top_n",      type = "integer", default = 5L),
  make_option("--min_pct",    type = "double",  default = 0.25),
  make_option("--logfc",      type = "double",  default = 0.25)
)))

res_col <- glue("clust_res{opt$resolution}")

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "RNA"

# Ensure a single joined data layer (required for v5 FindAllMarkers)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Force sequential future plan — parallel export of the test closure
# (~13 GB) exceeds future.globals.maxSize and silently skips all tests
future::plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)

Idents(obj) <- res_col
message(glue("[markers] FindAllMarkers on {nlevels(Idents(obj))} clusters ..."))

markers <- FindAllMarkers(obj, layer = "data", only.pos = TRUE,
                          min.pct = opt$min_pct, logfc.threshold = opt$logfc,
                          verbose = FALSE)
if (nrow(markers) == 0) stop("[markers] empty result — check layers and Seurat/SeuratObject versions")

write_csv(markers, "05_markers_all.csv")
message(glue("[markers] {nrow(markers)} rows written → 05_markers_all.csv"))

top <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = opt$top_n, with_ties = FALSE) |>
  ungroup()
write_csv(top, "05_markers_top.csv")

# Dotplot of deduplicated top markers
feats <- distinct(top, gene) |> pull(gene)
p <- DotPlot(obj, features = feats) + RotatedAxis() +
  labs(title = glue("Top {opt$top_n} markers / cluster ({res_col})")) +
  theme(axis.text.x = element_text(size = 7))
ggsave("05_dotplot_markers.pdf", p,
       width = max(12, 0.18 * length(feats)), height = 6)

message("[markers] done → 05_markers_top.csv + dotplot")
