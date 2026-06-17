#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — Step 03: VISUALIZE_QC
# -----------------------------------------------------------------------------
# Compute the remaining QC metrics (ribosomal %, transcriptome complexity) and
# render the per-sample QC plots used to DECIDE filtering thresholds. This step
# does NOT filter — it visualizes so you can eyeball each sample and, if a
# distribution differs from the defaults, set a qc_override in samples.yml
# before step 04 actually cuts.
#
# PLACE IN PIPELINE:
#   02 CREATE_SEURAT_OBJECT -> [03 VISUALIZE_QC] -> 04 FILTER_CELLS -> ...
#
# INPUT  : 02_seurat_unfiltered.rds (already has nFeature/nCount + percent.mt)
# OUTPUT : 03_seurat_with_qc.rds  (+ percent.ribo, log10GenesPerUMI)
#          03_qc_plots.pdf        (violins + scatters, current thresholds drawn)
# ENV    : paa-dev.yml  (Seurat 5.x, tidyverse, patchwork, fs, glue, optparse)
#
# The current per-sample thresholds are passed in only to be DRAWN as reference
# lines and to report how many cells they would remove — they are NOT applied.
# =============================================================================

# ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(tidyverse)   # ggplot2 / dplyr / readr
  library(patchwork)   # compose the multi-panel QC figure
  library(fs)
  library(glue)
})

# ---- 2. Command-line arguments ---------------------------------------------
# Thresholds are OPTIONAL: if omitted, plots render without reference lines.
option_list <- list(
  make_option("--in_rds",      type = "character", help = "02_seurat_unfiltered.rds"),
  make_option("--sample_id",   type = "character", help = "Sample id (plot titles / report)"),
  make_option("--ribo_pattern", type = "character", default = "^Rp[sl]",
              help = "Ribosomal gene-symbol pattern [default %default] (mouse)"),
  make_option("--min_features", type = "double", default = NA,
              help = "Lower nFeature threshold to DRAW (reference line only)"),
  make_option("--max_features", type = "double", default = NA,
              help = "Upper nFeature threshold to DRAW (reference line only)"),
  make_option("--max_percent_mt", type = "double", default = NA,
              help = "Max percent.mt to DRAW (reference line only)"),
  make_option("--out_rds",     type = "character", help = "Output .rds path"),
  make_option("--out_pdf",     type = "character", help = "Output QC plots .pdf path")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("in_rds", "sample_id", "out_rds", "out_pdf")
missing  <- keep(required, \(k) is.null(opt[[k]]))
if (length(missing) > 0)
  stop(glue("Missing required argument(s): {str_c(missing, collapse = ', ')}"))
if (!file_exists(opt$in_rds)) stop(glue("Input not found: {opt$in_rds}"))

# ---- 3. Load -----------------------------------------------------------------
message(glue("[03] Loading {opt$in_rds}"))
obj <- readRDS(opt$in_rds)

# ---- 4. Compute the remaining QC metrics ------------------------------------
# percent.mt already exists from step 02; add ribosomal % and complexity.
obj$percent.ribo <- PercentageFeatureSet(obj, pattern = opt$ribo_pattern)
# log10GenesPerUMI = transcriptome "complexity": low values flag low-diversity
# cells (e.g. dying cells, red blood cells) even when counts look adequate.
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

meta <- as_tibble(obj@meta.data)

# ---- 5. Plot helpers --------------------------------------------------------
# One violin per metric, with optional dashed reference line(s) at thresholds.
qc_violin <- function(df, metric, hlines = NULL) {
  p <- ggplot(df, aes(x = sample_id, y = .data[[metric]])) +
    geom_violin(fill = "grey85", colour = "grey40") +
    geom_jitter(width = 0.2, size = 0.1, alpha = 0.2) +
    labs(x = NULL, y = metric) +
    theme_bw()
  hlines <- hlines[!is.na(hlines)]
  if (length(hlines) > 0)
    p <- p + geom_hline(yintercept = hlines, linetype = "dashed", colour = "firebrick")
  p
}
# Bivariate QC scatter coloured by a third metric (the doublet/quality view).
qc_scatter <- function(df, x, y, colour = "percent.mt") {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_colour_viridis_c() +
    theme_bw()
}

# ---- 6. Assemble and save the QC figure -------------------------------------
violins <- qc_violin(meta, "nFeature_RNA", c(opt$min_features, opt$max_features)) |
           qc_violin(meta, "nCount_RNA") |
           qc_violin(meta, "percent.mt", opt$max_percent_mt) |
           qc_violin(meta, "percent.ribo")

scatters <- qc_scatter(meta, "nCount_RNA", "nFeature_RNA") |
            qc_scatter(meta, "nCount_RNA", "percent.mt")

qc_fig <- (violins / scatters) +
  plot_annotation(
    title    = glue("QC — {opt$sample_id}"),
    subtitle = glue("{ncol(obj)} cells  |  dashed red = current filter thresholds (step 04)")
  )

dir_create(path_dir(opt$out_pdf))
ggsave(opt$out_pdf, qc_fig, width = 14, height = 9)
message(glue("[03] Wrote {opt$out_pdf}"))

# ---- 7. Report distributions + threshold-impact preview ---------------------
message(glue(
  "[03] QC medians for '{opt$sample_id}':\n",
  "       nFeature_RNA : {round(median(meta$nFeature_RNA))}\n",
  "       nCount_RNA   : {round(median(meta$nCount_RNA))}\n",
  "       percent.mt   : {round(median(meta$percent.mt), 2)} %\n",
  "       percent.ribo : {round(median(meta$percent.ribo), 2)} %"
))
# If thresholds were supplied, report how many cells they WOULD remove. This is
# the decision-support number: a high % is the cue to revisit the defaults.
if (!anyNA(c(opt$min_features, opt$max_features, opt$max_percent_mt))) {
  fail <- meta$nFeature_RNA < opt$min_features |
          meta$nFeature_RNA > opt$max_features |
          meta$percent.mt   > opt$max_percent_mt
  message(glue(
    "[03] At current thresholds (nFeature {opt$min_features}-{opt$max_features}, ",
    "mt <= {opt$max_percent_mt}%): {sum(fail)}/{length(fail)} cells would be removed ",
    "({round(100 * mean(fail), 1)}%)."
  ))
}

# ---- 8. Save the QC-annotated object ----------------------------------------
saveRDS(obj, opt$out_rds)
message(glue("[03] Wrote {opt$out_rds}"))