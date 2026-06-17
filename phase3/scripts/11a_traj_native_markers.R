#!/usr/bin/env Rscript
# 11a_traj_native_markers.R
# The 09 CSV/dotplot artifacts sitting in endo_traj/ are STALE endo+endocardium
# copies (summary totals 1552 cells with an Endocardium lineage; clusters 1/4/7 are
# endocardial by marker). Before annotating the native uncorrected clusters, this
# script (a) verifies the actual endo_traj object's composition and (b) regenerates
# markers + timepoint composition + a dotplot FROM the object, so the annotation is
# built on files that provably match the cells the trajectory runs on.
#
# If the verification reports ~1552 cells / endocardial markers present, the RDS
# itself still contains endocardium and 09 (--integrate FALSE) must be re-run on
# endothelium-only lineages before proceeding. If it reports ~1049 / endocardial
# markers absent, the uploaded files were merely stale and these outputs are correct.
#
# Input : results/phase2a/endo_traj/09_endo_subset.rds
# Output: 11a_native_markers_top.csv | 11a_native_summary.csv | 11a_native_dotplot.pdf

suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(glue); library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",  type = "character",
              default = "results/phase2a/endo_traj/09_endo_subset.rds"),
  make_option("--outdir", type = "character", default = "results/phase2a/endo_traj")
)))
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
  write_csv(comp, out("11a_native_summary.csv")); print(comp)
}

# -- (b) Markers per native cluster ---------------------------------------------
mk  <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top <- mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 8) |> ungroup()
write_csv(top, out("11a_native_markers_top.csv"))

# -- (b) Dotplot on the bed-marker panel ----------------------------------------
panel <- c("Etv2","Tal1","Lmo2","Fli1","Kdr","Pecam1","Cdh5","Cldn5","Emcn",
           "Gja5","Gja4","Efnb2","Dll4","Hey1","Sox17","Nr2f2","Nrp2","Aplnr","Ephb4",
           "Npr3","Nfatc1","Gata5","Klf2","Bmx","Prox1","Lyve1","Flt4",
           "Cxcr4","Itga5","Itgb1","Fn1")
panel <- panel[panel %in% rownames(obj)]
dp <- DotPlot(obj, features = panel) + RotatedAxis() +
  ggtitle("endo_traj native clusters — bed marker panel")
ggsave(out("11a_native_dotplot.pdf"), dp, width = 14, height = 5)

message(glue("[done] -> {out('11a_native_markers_top.csv')}"))