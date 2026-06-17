#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 04: FILTER_CELLS
# -----------------------------------------------------------------------------
# Apply the per-sample QC thresholds decided after inspecting step 03's plots,
# and report exactly how many cells each criterion removes (transparent, so the
# filtering is auditable rather than a silent cell drop).
#
# PLACE IN PIPELINE:
#   03 VISUALIZE_QC -> [04 FILTER_CELLS] -> 05 DETECT_DOUBLETS -> ...
#
# INPUT  : 03_seurat_with_qc.rds
# OUTPUT : 04_seurat_filtered.rds   (feeds DETECT_DOUBLETS at step 05)
# ENV    : paa-dev.yml  (Seurat 5.x, tidyverse, fs, glue, optparse)
#
# SCOPE NOTE: this filters on EXACTLY the three parameterised thresholds
# (nFeature lower/upper, percent.mt upper) — nothing else, so the cut is fully
# described by samples.yml. The max_features ceiling here is a COARSE upper
# bound (extreme outliers only); proper doublet identification is the job of
# step 05 (scDblFinder), so max_features can be set generously upstream.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(tidyverse)   # purrr::keep / stringr::str_c in the arg check
  library(fs)
  library(glue)
})

# ---- 2. Command-line arguments ---------------------------------------------
# Flags mirror exactly what phase1/main.nf passes in FILTER_CELLS.
option_list <- list(
  make_option("--in_rds",         type = "character", help = "03_seurat_with_qc.rds"),
  make_option("--min_features",   type = "double", help = "Minimum nFeature_RNA to KEEP"),
  make_option("--max_features",   type = "double", help = "Maximum nFeature_RNA to KEEP (coarse outlier cap)"),
  make_option("--max_percent_mt", type = "double", help = "Maximum percent.mt to KEEP"),
  make_option("--out",            type = "character", help = "Output 04_seurat_filtered.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("in_rds", "min_features", "max_features", "max_percent_mt", "out")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))
if (!file_exists(opt$in_rds)) stop(glue("Input not found: {opt$in_rds}"))

# ---- 3. Load ----------------------------------------------------------------
message(glue("[04] Loading {opt$in_rds}"))
obj <- readRDS(opt$in_rds)
md  <- obj@meta.data
n0  <- nrow(md)
sid <- md$sample_id[1]

# ---- 4. Build the keep-mask from metadata (robust; avoids subset() NSE) -----
# Per-criterion failure masks so we can report each one (they overlap; a cell
# can fail more than one, so the per-criterion counts need not sum to the total).
fail_low  <- md$nFeature_RNA <  opt$min_features
fail_high <- md$nFeature_RNA >  opt$max_features
fail_mt   <- md$percent.mt   >  opt$max_percent_mt
keep_mask <- !(fail_low | fail_high | fail_mt)

if (sum(keep_mask) == 0)
  stop(glue(
    "All {n0} cells fail the filter for '{sid}' — thresholds are almost certainly ",
    "wrong (nFeature {opt$min_features}-{opt$max_features}, mt <= {opt$max_percent_mt}%)."
  ))

# ---- 5. Subset to retained cells --------------------------------------------
obj <- subset(obj, cells = colnames(obj)[keep_mask])
n1  <- ncol(obj)

# ---- 6. Transparent removal report (captured in the Nextflow log) -----------
message(glue(
  "[04] Filtering '{sid}':\n",
  "       fail nFeature < {opt$min_features}     : {sum(fail_low)}\n",
  "       fail nFeature > {opt$max_features}     : {sum(fail_high)}\n",
  "       fail percent.mt > {opt$max_percent_mt}%   : {sum(fail_mt)}\n",
  "       cells: {n0} -> {n1}  ({round(100 * n1 / n0, 1)}% retained, {n0 - n1} removed)"
))
# Sanity guard: a very low retention usually means a misconfigured threshold for
# this particular sample — the cue to set a qc_override in samples.yml.
if (n1 / n0 < 0.10)
  warning(glue("Only {round(100 * n1 / n0, 1)}% of cells retained for '{sid}' — ",
               "verify this sample's thresholds (consider a qc_override)."))

# ---- 7. Save ----------------------------------------------------------------
saveRDS(obj, opt$out)
message(glue("[04] Wrote {opt$out}"))