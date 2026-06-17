#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 02: CREATE_SEURAT_OBJECT
# -----------------------------------------------------------------------------
# Build a per-sample Seurat object from the canonical 10x triplet that
# scripts/download_geo.sh placed in raw_data/<id>/, attach sample metadata,
# and compute the mitochondrial-% QC metric. NO cell filtering happens here —
# that is deferred to step 04 (single responsibility: 02 builds, 04 filters).
#
# PLACE IN PIPELINE:
#   (download_geo.sh) -> [02 CREATE_SEURAT_OBJECT] -> 03 VISUALIZE_QC -> ...
#
# INPUT  : raw_data/<id>/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}
# OUTPUT : 02_seurat_unfiltered.rds  (raw object + metadata + percent.mt)
# ENV    : paa-dev.yml  (Seurat 5.x, tidyverse, fs, glue, optparse)
#
# FEATURE RECONCILIATION NOTE:
#   A per-sample script cannot intersect gene sets across samples — that is a
#   merge-time job (Phase 3a). What we DO here to make that intersection safe
#   later: preserve the full features table (Ensembl ID + symbol + type) in
#   obj@misc$feature_metadata, so the merge can reconcile on STABLE Ensembl IDs
#   rather than on symbols (which can collide or drift). All 7 samples are
#   mm10 / Cell Ranger 3.1.0, so the gene lists should already be near-identical.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(tidyverse)   # readr / purrr / stringr used below
  library(fs)          # modern, consistent path handling
  library(glue)        # readable string interpolation in messages
})

# ---- 2. Command-line arguments ---------------------------------------------
# Flags mirror exactly what phase1/main.nf passes in CREATE_SEURAT_OBJECT.
option_list <- list(
  make_option("--raw_dir",    type = "character", help = "Dir with the 10x triplet (raw_data/<id>/)"),
  make_option("--sample_id",  type = "character", help = "Sample id, e.g. mesp1_ko_e95_rep1"),
  make_option("--genotype",   type = "character", help = "WT or KO"),
  make_option("--timepoint",  type = "character", help = "E80 / E825 / E95 / E105"),
  make_option("--series_id",  type = "character", help = "GEO SubSeries id (batch key)"),
  make_option("--mt_pattern", type = "character", default = "^mt-",
              help = "Mitochondrial gene-symbol prefix [default %default] (mouse)"),
  make_option("--min_cells",  type = "integer", default = 3L,
              help = "Drop genes detected in fewer than this many cells [default %default]"),
  make_option("--out",        type = "character", help = "Output .rds path")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Fail early and clearly if a required argument is missing.
required <- c("raw_dir", "sample_id", "genotype", "timepoint", "series_id", "out")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))

# ---- 3. Validate inputs -----------------------------------------------------
# Read10X expects this exact canonical triplet (download_geo.sh guarantees it).
need   <- path(opt$raw_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
absent <- need[!file_exists(need)]
if (length(absent) > 0)
  stop(glue(
    "Missing 10x file(s) in {opt$raw_dir}:\n  ",
    "{str_c(path_file(absent), collapse = '\n  ')}\n",
    "(Did you run scripts/download_geo.sh first?)"
  ))

# ---- 4. Read counts and build the object ------------------------------------
message(glue("[02] Reading 10x matrix for '{opt$sample_id}' from {opt$raw_dir}"))
counts <- Read10X(data.dir = opt$raw_dir)            # gene.column = 2 (symbols) by default
# Read10X returns a list when several feature types exist (e.g. Antibody Capture);
# these samples are gene expression only, but handle the list case defensively.
if (is.list(counts)) {
  counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
}

# min.features = 0 here ON PURPOSE: all cell-level QC lives in step 04. min.cells
# only drops genes seen in <min_cells cells (shrinks the object; not a cell filter).
obj <- CreateSeuratObject(
  counts       = counts,
  project      = opt$sample_id,
  min.cells    = opt$min_cells,
  min.features = 0
)

# ---- 5. Attach sample metadata ----------------------------------------------
# These travel with every cell and drive grouping/blocking downstream
# (genotype = the contrast; series_id = the batch key for integration).
obj$sample_id <- opt$sample_id
obj$genotype  <- opt$genotype
obj$timepoint <- opt$timepoint
obj$series_id <- opt$series_id

# ---- 6. Mitochondrial percentage (the key QC metric used by step 04) --------
obj$percent.mt <- PercentageFeatureSet(obj, pattern = opt$mt_pattern)
# Guard: if NOTHING matched the mt pattern, percent.mt is all zero — usually a
# sign the rownames are Ensembl IDs, not symbols, or the wrong species pattern.
if (all(obj$percent.mt == 0))
  warning(glue("percent.mt is 0 for ALL cells — check --mt_pattern ('{opt$mt_pattern}') ",
               "and that rownames are gene symbols."))

# ---- 7. Preserve the features table for merge-time reconciliation -----------
# Stored in misc (not as rownames) so later phases can intersect on Ensembl IDs.
feat <- tryCatch(
  read_tsv(path(opt$raw_dir, "features.tsv.gz"),
           col_names = FALSE, show_col_types = FALSE),
  error = function(e) NULL
)
if (!is.null(feat)) {
  names(feat) <- c("gene_id", "gene_symbol", "feature_type")[seq_len(ncol(feat))]
  obj@misc$feature_metadata <- feat
  style <- if (ncol(feat) >= 3) "CR3-style features.tsv" else "CR2-style genes.tsv"
  message(glue("[02] Features table: {nrow(feat)} genes, {ncol(feat)} columns ({style})"))
}

# ---- 8. Report (captured in the Nextflow process log) -----------------------
message(glue(
  "[02] Built object for '{opt$sample_id}':\n",
  "       cells              : {ncol(obj)}\n",
  "       genes              : {nrow(obj)}\n",
  "       median genes/cell  : {round(median(obj$nFeature_RNA))}\n",
  "       median counts/cell : {round(median(obj$nCount_RNA))}\n",
  "       median percent.mt  : {round(median(obj$percent.mt), 2)} %"
))

# ---- 9. Save ----------------------------------------------------------------
saveRDS(obj, opt$out)
message(glue("[02] Wrote {opt$out}"))