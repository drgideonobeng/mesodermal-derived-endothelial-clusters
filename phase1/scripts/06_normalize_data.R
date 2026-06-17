#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 06: NORMALIZE_DATA  (SCTransform)
# -----------------------------------------------------------------------------
# Variance-stabilizing normalization of the doublet-free singlets, regressing
# out mitochondrial fraction, and selection of variable features for PCA.
#
# PLACE IN PIPELINE:
#   05 DETECT_DOUBLETS -> [06 NORMALIZE_DATA] -> 07 DIM_REDUCTION_AND_CLUSTER
#
# INPUT  : 05_seurat_singlets.rds
# OUTPUT : 06_seurat_normalized.rds   (adds an 'SCT' assay, set as default)
# ENV    : paa-dev.yml  (Seurat 5.x, glmGamPoi, tidyverse, fs, glue, optparse)
#
# NOTES:
#   * vst.flavor = "v2" is the current SCTransform (uses the glmGamPoi backend
#     in the env — faster and the recommended default in Seurat 5).
#   * vars.to.regress = "percent.mt" removes residual mitochondrial signal, per
#     the mesp1/tbx1 design guide.
#   * SCTransform is per-sample here ON PURPOSE — its noise model is dataset-
#     specific, so it must NOT be fit across the merged multi-series object.
#     Cross-sample integration (Phase 3a) starts from these per-sample objects.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(tidyverse)
  library(fs)
  library(glue)
})

# ---- 2. Command-line arguments ---------------------------------------------
option_list <- list(
  make_option("--in_rds", type = "character", help = "05_seurat_singlets.rds"),
  make_option("--n_hvgs", type = "integer", default = 2000L,
              help = "Number of variable features to select [default %default]"),
  make_option("--seed",   type = "integer", default = 42L,
              help = "RNG seed for SCTransform [default %default]"),
  make_option("--out",    type = "character", help = "Output 06_seurat_normalized.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("in_rds", "out")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))
if (!file_exists(opt$in_rds)) stop(glue("Input not found: {opt$in_rds}"))

# ---- 3. Load ----------------------------------------------------------------
message(glue("[06] Loading {opt$in_rds}"))
obj <- readRDS(opt$in_rds)
sid <- obj$sample_id[1]

# ---- 4. SCTransform ---------------------------------------------------------
# SCTransform uses the `future` package for internal parallelism and serializes
# a large function closure (~500 MB) to each worker. Raise future's per-worker
# global-size cap so it does not reject this as a potential mistake.
# This is expected SCTransform behaviour — not a memory leak.
options(future.globals.maxSize = 8000 * 1024^2)   # 8 GB

set.seed(opt$seed)
obj <- SCTransform(
  obj,
  assay               = "RNA",
  vars.to.regress     = "percent.mt",
  variable.features.n = opt$n_hvgs,
  vst.flavor          = "v2",
  seed.use            = opt$seed,
  verbose             = FALSE
)

# ---- 5. Report --------------------------------------------------------------
message(glue(
  "[06] SCTransform complete for '{sid}':\n",
  "       default assay     : {DefaultAssay(obj)}\n",
  "       variable features : {length(VariableFeatures(obj))}\n",
  "       cells x genes     : {ncol(obj)} x {nrow(obj)}"
))

# ---- 6. Save ----------------------------------------------------------------
saveRDS(obj, opt$out)
message(glue("[06] Wrote {opt$out}"))