#!/usr/bin/env Rscript
# 01_normalize.R
# Merge four pre-QC Seurat objects, log-normalise, select HVGs.
# Input : samplesheet CSV (sample_id | path | timepoint)
# Output: 01_merged.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--samplesheet", type = "character"),
  make_option("--n_hvgs",      type = "integer", default = 2000L)
)))

tp_order <- c("E80", "E825", "E95", "E105")

# ── Load samples ──────────────────────────────────────────────────────────────
ss   <- read_csv(opt$samplesheet, show_col_types = FALSE)
objs <- map2(ss$path, ss$sample_id, function(path, id) {
  obj             <- readRDS(path)
  obj$sample_id   <- id
  obj$timepoint   <- factor(ss$timepoint[ss$sample_id == id], levels = tp_order)
  obj
})
message(glue("[normalize] loaded {nrow(ss)} samples"))

# ── Merge ─────────────────────────────────────────────────────────────────────
merged <- merge(objs[[1]], objs[-1], add.cell.ids = ss$sample_id) |>
  JoinLayers()
message(glue("[normalize] {ncol(merged)} cells after merge"))

# ── Normalise + HVG selection ─────────────────────────────────────────────────
merged <- merged |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = opt$n_hvgs, verbose = FALSE)

message(glue("[normalize] {opt$n_hvgs} HVGs selected"))
saveRDS(merged, "01_merged.rds")
message("[normalize] done → 01_merged.rds")
