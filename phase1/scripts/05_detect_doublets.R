#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 05: DETECT_DOUBLETS  (scDblFinder, per sample)
# -----------------------------------------------------------------------------
# Identify and remove droplet doublets within a single sample, BEFORE
# normalization and integration. Runs on the QC-filtered counts from step 04.
#
# PLACE IN PIPELINE:
#   04 FILTER_CELLS -> [05 DETECT_DOUBLETS] -> 06 NORMALIZE_DATA -> ...
#   The singlets object (05) is the hand-off to Phase 3a integration — clean
#   singlets must be what gets integrated.
#
# INPUT  : 04_seurat_filtered.rds
# OUTPUT : 05_seurat_singlets.rds   (doublets removed; class+score kept in meta)
#          05_doublet_summary.csv   (per-sample rate)
#          05_doublet_plots.pdf     (score distribution + count-space diagnostic)
# ENV    : paa-dev.yml  (Seurat 5.x, scDblFinder, SingleCellExperiment,
#                        tidyverse, patchwork, fs, glue, optparse)
#
# WHY RANDOM DOUBLET GENERATION (use_clusters = FALSE by default):
#   scDblFinder builds artificial doublets to train its classifier. Cluster-
#   based generation assumes discrete cell types; on data with a continuous
#   trajectory (our SHF -> angioblast -> arterial-EC axis) it is more likely to
#   score genuine INTERMEDIATE cells as doublets — deleting the very biology we
#   want. Random generation is the trajectory-safe choice here. Flip to TRUE
#   only for samples with clearly discrete structure and no continuum of interest.
#
# HONEST CAVEAT: heterotypic doublets and real transitional cells look alike
#   (both co-express two programs). This step FLAGS then removes, and keeps the
#   score so you can judge the call: inspect the score histogram — a clean
#   bimodal split is safe; a smooth continuum is a sign intermediates may be
#   getting scored as doublets, in which case loosen / revisit.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(optparse)
  library(tidyverse)
  library(patchwork)
  library(fs)
  library(glue)
})

# ---- 2. Command-line arguments ---------------------------------------------
option_list <- list(
  make_option("--in_rds",      type = "character", help = "04_seurat_filtered.rds"),
  make_option("--sample_id",   type = "character", help = "Sample id (report / plots)"),
  make_option("--seed",        type = "integer", default = 42L,
              help = "RNG seed — scDblFinder simulates doublets stochastically [default %default]"),
  make_option("--use_clusters", type = "logical", default = FALSE,
              help = "Cluster-based (TRUE) vs random (FALSE) doublet generation [default %default]"),
  make_option("--out_rds",     type = "character", help = "Output 05_seurat_singlets.rds"),
  make_option("--out_summary", type = "character", help = "Output 05_doublet_summary.csv"),
  make_option("--out_pdf",     type = "character", help = "Output 05_doublet_plots.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("in_rds", "sample_id", "out_rds", "out_summary", "out_pdf")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))
if (!file_exists(opt$in_rds)) stop(glue("Input not found: {opt$in_rds}"))

# ---- 3. Load ----------------------------------------------------------------
message(glue("[05] Loading {opt$in_rds}"))
obj <- readRDS(opt$in_rds)

# ---- 4. Run scDblFinder on the raw counts -----------------------------------
# Pass the count matrix directly (avoids Seurat<->SCE object-conversion quirks);
# scDblFinder returns an SCE carrying the per-cell class and score.
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
set.seed(opt$seed)
sce <- scDblFinder(counts, clusters = opt$use_clusters)

# Cell-order guard: results must align 1:1 with the Seurat object before we
# copy them back (a mismatch here would mislabel cells — fail loudly instead).
stopifnot("scDblFinder cell order differs from object" =
            identical(colnames(sce), colnames(obj)))
obj$scDblFinder.class <- sce$scDblFinder.class
obj$scDblFinder.score <- sce$scDblFinder.score

# ---- 5. Summary table -------------------------------------------------------
n_cells   <- ncol(obj)
n_doublet <- sum(obj$scDblFinder.class == "doublet")
summary_tbl <- tibble(
  sample_id    = opt$sample_id,
  n_cells      = n_cells,
  n_doublet    = n_doublet,
  n_singlet    = n_cells - n_doublet,
  doublet_rate = round(n_doublet / n_cells, 4),
  generation   = if (opt$use_clusters) "cluster-based" else "random"
)
dir_create(path_dir(opt$out_summary))
write_csv(summary_tbl, opt$out_summary)

# ---- 6. Diagnostic plots ----------------------------------------------------
md <- as_tibble(obj@meta.data)
# (a) Score distribution by called class — judge whether the split is clean.
p_score <- ggplot(md, aes(x = scDblFinder.score, fill = scDblFinder.class)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  labs(title = "Doublet score distribution", x = "scDblFinder.score", y = "cells", fill = NULL) +
  theme_bw()
# (b) Count space — true heterotypic doublets tend to sit at high counts/genes.
p_counts <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA, colour = scDblFinder.class)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "Doublets in count space", colour = NULL) +
  theme_bw()

qc_fig <- (p_score | p_counts) +
  plot_annotation(
    title    = glue("Doublet detection — {opt$sample_id}"),
    subtitle = glue("{n_doublet}/{n_cells} called doublet ",
                    "({round(100 * n_doublet / n_cells, 1)}%)  |  ",
                    "{summary_tbl$generation} generation")
  )
dir_create(path_dir(opt$out_pdf))
ggsave(opt$out_pdf, qc_fig, width = 12, height = 5)
message(glue("[05] Wrote {opt$out_pdf} and {opt$out_summary}"))

# ---- 7. Remove doublets -> singlets object ----------------------------------
singlets <- subset(obj, cells = colnames(obj)[obj$scDblFinder.class == "singlet"])
message(glue(
  "[05] '{opt$sample_id}': {n_doublet}/{n_cells} doublets removed ",
  "({round(100 * n_doublet / n_cells, 1)}%) -> {ncol(singlets)} singlets. ",
  "Inspect {path_file(opt$out_pdf)}: clean bimodal score = safe; smooth ",
  "continuum = intermediates may be over-called."
))

# ---- 8. Save ----------------------------------------------------------------
saveRDS(singlets, opt$out_rds)
message(glue("[05] Wrote {opt$out_rds}"))
