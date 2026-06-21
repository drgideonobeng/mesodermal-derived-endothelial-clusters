#!/usr/bin/env Rscript
# 03a_traj_native_markers.R
# Verify the uncorrected object's composition and regenerate native-cluster markers
# directly FROM that object, so the native (uncorrected) cluster annotation is built
# on files that provably match the cells the trajectory runs on — not on stale CSVs.
#
# History note: the original endo_traj/ CSV+dotplot artifacts were STALE
# endo+endocardium copies (summary totalled 1552 cells with an Endocardium lineage;
# clusters 1/4/7 were endocardial by marker). This script (a) verifies the actual
# uncorrected object's composition and (b) regenerates markers + timepoint
# composition + a dotplot from the object.
#
# If verification reports ~1552 cells / endocardial markers present, the RDS still
# contains endocardium and 01 (--integrate FALSE) must be re-run on endothelium-only
# lineages before proceeding. If it reports ~1049 / endocardial markers absent, the
# object is endothelium-only and these outputs are correct.
#
# Input : results/phase3/endothelium/uncorrected/01_endo_subset.rds
# Output: 03a_native_markers_top.csv | 03a_native_summary.csv | 03a_native_dotplot.pdf

suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(glue); library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage", type = "character", default = "endothelium",
              help = "phase-3 namespace under results/phase3/<lineage>/"),
  make_option("--input",   type = "character", default = NULL,
              help = "uncorrected subset; default results/phase3/<lineage>/uncorrected/01_endo_subset.rds"),
  make_option("--outdir",  type = "character", default = NULL,
              help = "default results/phase3/<lineage>/uncorrected")
)))

# ── Resolve paths from --lineage (override individually if needed) ─────────────
base <- file.path("results/phase3", opt$lineage)
if (is.null(opt$input))  opt$input  <- file.path(base, "uncorrected", "01_endo_subset.rds")
if (is.null(opt$outdir)) opt$outdir <- file.path(base, "uncorrected")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)   # Seurat v5 split-layer guard

# -- (a) VERIFY composition -----------------------------------------------------
message(glue("[verify] cells: {ncol(obj)}   (expect 1049 = endothelium-only)"))
lin_col <- intersect(c("lineage","prior_lineage","lineage_07","celltype"),
                     colnames(obj@meta.data))[1]
if (!is.na(lin_col)) { message("[verify] prior-lineage make-up:"); print(table(obj@meta.data[[lin_col]])) }

endo_mk <- intersect(c("Npr3","Nfatc1","Gata5","Bmx"), rownames(obj))
if (length(endo_mk)) {
  pct <- colMeans(FetchData(obj, endo_mk) > 0) * 100
  message("[verify] %cells+ for endocardial markers (high => endocardium still present):")
  print(round(pct, 1))
}

# -- Identify native clustering + timepoint columns -----------------------------
clu <- intersect(c("subcluster","seurat_clusters"), colnames(obj@meta.data))[1]
tp  <- intersect(c("timepoint","stage","orig.ident"), colnames(obj@meta.data))[1]
stopifnot(!is.na(clu))
Idents(obj) <- clu
message(glue("[clusters] native '{clu}': {paste(levels(Idents(obj)), collapse = ', ')}"))

# -- (b) Timepoint composition per native cluster -------------------------------
if (!is.na(tp)) {
  comp <- obj@meta.data |> as_tibble() |>
    count(subcluster = .data[[clu]], timepoint = .data[[tp]], name = "n")
  write_csv(comp, out("03a_native_summary.csv")); print(comp)
}

# -- (b) Markers per native cluster ---------------------------------------------
mk  <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top <- mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 8) |> ungroup()
write_csv(top, out("03a_native_markers_top.csv"))

# -- (b) Dotplot on the bed-marker panel ----------------------------------------
panel <- c("Etv2","Tal1","Lmo2","Fli1","Kdr","Pecam1","Cdh5","Cldn5","Emcn",
           "Gja5","Gja4","Efnb2","Dll4","Hey1","Sox17","Nr2f2","Nrp2","Aplnr","Ephb4",
           "Npr3","Nfatc1","Gata5","Klf2","Bmx","Prox1","Lyve1","Flt4",
           "Cxcr4","Itga5","Itgb1","Fn1")
panel <- panel[panel %in% rownames(obj)]
dp <- DotPlot(obj, features = panel) + RotatedAxis() +
  ggtitle("uncorrected native clusters — bed marker panel")
ggsave(out("03a_native_dotplot.pdf"), dp, width = 14, height = 5)

message(glue("[done] -> {out('03a_native_markers_top.csv')}"))
