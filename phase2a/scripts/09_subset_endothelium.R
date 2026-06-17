#!/usr/bin/env Rscript
# 09_subset_endothelium.R
# Subset the endothelial compartment BY LINEAGE LABEL (not by single-gene gating)
# and RE-DERIVE structure within it: re-HVG, re-integrate, re-cluster. The point
# is to resolve the endothelial *beds* (arterial / venous / endocardial /
# hemogenic / PAA) that sit below the res-0.4 whole-embryo clustering, so the
# subsequent trajectory step has named endpoints and a real root to work with.
#
# ── WHY SUBSET BY LABEL, NOT BY GENE ──────────────────────────────────────────
# Cell identity is multivariate; the lineage labels are a robust boundary. A
# single-gene gate (e.g. "Cldn5 > 0") is corrupted by ambient transcript and
# dropout in both directions. We subset on `lineage`, the frozen identity.
#
# ── THE INTEGRATION FORK (read before running) ────────────────────────────────
# `timepoint` here is BOTH the batch variable AND the biological time axis.
#   --integrate TRUE  (default): fastMNN mixes timepoints out. Use this to
#       re-CLUSTER and find beds — you want bed identity, not stage, to drive
#       the clusters.
#   --integrate FALSE: uncorrected PCA, timepoints preserved. Use this to build
#       the embedding you'll hand to TRAJECTORY — fastMNN can flatten the very
#       developmental gradient pseudotime needs to trace.
# Run both; compare the by-timepoint UMAP. If the integrated version fully
# homogenises timepoints, that's over-correction for trajectory purposes.
#
# ── OTHER NOTES ───────────────────────────────────────────────────────────────
# • HVGs are RE-SELECTED on the subset. Whole-embryo HVGs encode between-lineage
#   variation that is constant inside the endothelial compartment; reusing them
#   would just recover Endothelium-vs-Endocardium and nothing finer.
# • Small subset (Endocardium esp. is a small cluster). Small n over-fragments;
#   validate every subcluster against markers, do not trust the count. Use
#   --run_sweep to pick resolution by stability rather than guessing.
# • Cell cycle was scored but not regressed upstream. A cycling sub-state can
#   split out and masquerade as a bed. --exclude_cc drops Tirosh CC genes from
#   the HVG list to suppress that (affects the embedding; off by default).
# • fastMNN needs enough cells per batch. Timepoints with < --min_cells_batch
#   cells are dropped from the integration (reported), since they can't be
#   reliably corrected.
#
# Input : results/phase2a/objects/07_annotated.rds  (frozen `lineage`)
# Output: 09_endo_subset.rds | 09_endo_umap_overlays.pdf |
#         09_endo_bed_dotplot.pdf | 09_endo_subcluster_markers.csv |
#         09_endo_subcluster_top.csv | 09_endo_subcluster_summary.csv
#         (+ 09_endo_clustree.pdf if --run_sweep)

suppressPackageStartupMessages({
  library(Seurat)
  library(batchelor)
  library(BiocParallel)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",           type = "character",
              default = "results/phase2a/objects/07_annotated.rds"),
  make_option("--lineages",        type = "character",
              default = "Endothelium,Endocardium",
              help = "comma-separated frozen lineage labels to subset"),
  make_option("--n_hvgs",          type = "integer", default = 2000L),
  make_option("--n_dims",          type = "integer", default = 30L),
  make_option("--resolution",      type = "double",  default = 0.6),
  make_option("--k",               type = "integer", default = 20L,
              help = "fastMNN nearest neighbours; lower for small subsets"),
  make_option("--min_cells_batch", type = "integer", default = 20L),
  make_option("--integrate",       type = "logical", default = TRUE,
              help = "TRUE: fastMNN (for bed clustering). FALSE: PCA (for trajectory)"),
  make_option("--exclude_cc",      action = "store_true", default = FALSE,
              help = "drop Tirosh cell-cycle genes from the HVG list"),
  make_option("--run_sweep",       action = "store_true", default = FALSE,
              help = "clustree resolution sweep (does not change locked resolution)"),
  make_option("--outdir",          type = "character", default = "results/phase2a/endo_subset",
              help = "directory for all outputs (created if absent)")
)))

set.seed(42)
tp_order <- c("E80", "E825", "E95", "E105")
tp_cols  <- setNames(viridisLite::viridis(4, option = "D", direction = -1), tp_order)
lineages <- str_split(opt$lineages, ",")[[1]] |> str_trim()

# All outputs go under --outdir (bare filenames otherwise land in the CWD,
# since this script emits unprefixed names for Nextflow publishDir routing).
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

# ── Subset by lineage label ───────────────────────────────────────────────────
obj <- readRDS(opt$input)
stopifnot("lineage" %in% colnames(obj@meta.data))
Idents(obj) <- "lineage"

absent <- setdiff(lineages, levels(droplevels(obj$lineage)))
if (length(absent) > 0)
  stop(glue("[subset] lineages not found in object: {paste(absent, collapse=', ')}"))

sub <- subset(obj, idents = lineages)
sub$lineage <- droplevels(sub$lineage)
message(glue("[subset] {ncol(sub)} cells across {paste(lineages, collapse=' + ')}"))
print(table(prior_lineage = sub$lineage, timepoint = sub$timepoint))

# ── Drop tiny batches (integration cannot correct them) ───────────────────────
tp_n   <- table(factor(as.character(sub$timepoint), levels = tp_order))
tp_keep <- names(tp_n)[tp_n >= opt$min_cells_batch]
tp_drop <- setdiff(names(tp_n)[tp_n > 0], tp_keep)
if (length(tp_drop) > 0 && opt$integrate) {
  message(glue("[subset] dropping under-populated timepoints from integration: ",
               "{paste(glue('{tp_drop} (n={tp_n[tp_drop]})'), collapse=', ')}"))
  sub <- subset(sub, subset = timepoint %in% tp_keep)
}
tp_present <- tp_order[tp_order %in% as.character(unique(sub$timepoint))]

# ── Re-normalise + RE-SELECT HVGs on the subset ───────────────────────────────
DefaultAssay(sub) <- "RNA"
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
sub <- sub |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = opt$n_hvgs, verbose = FALSE)

if (opt$exclude_cc) {
  cc <- c(str_to_title(cc.genes.updated.2019$s.genes),
          str_to_title(cc.genes.updated.2019$g2m.genes))
  vf <- setdiff(VariableFeatures(sub), cc)
  VariableFeatures(sub) <- vf
  message(glue("[subset] excluded {opt$n_hvgs - length(vf)} cell-cycle genes from HVGs"))
}
hvg <- VariableFeatures(sub)

# ── Embedding: fastMNN (bed clustering) OR PCA (trajectory) ────────────────────
if (opt$integrate && length(tp_present) >= 2) {
  message(glue("[subset] fastMNN over timepoints: {paste(tp_present, collapse=' -> ')} (k={opt$k})"))
  log_mat <- GetAssayData(sub, layer = "data")
  mnn_out <- fastMNN(log_mat, batch = factor(as.character(sub$timepoint), levels = tp_present),
                     subset.row = hvg, d = opt$n_dims, k = opt$k,
                     merge.order = tp_present, BPPARAM = SerialParam())
  emb <- reducedDim(mnn_out, "corrected"); rownames(emb) <- colnames(sub)
  sub[["mnn_sub"]] <- CreateDimReducObject(emb, key = "MNNSUB_", assay = "RNA")
  red <- "mnn_sub"
} else {
  if (opt$integrate) message("[subset] < 2 usable timepoints — falling back to PCA")
  message("[subset] uncorrected PCA embedding (timepoints preserved)")
  sub <- sub |> ScaleData(verbose = FALSE) |>
    RunPCA(npcs = opt$n_dims, reduction.name = "pca_sub", verbose = FALSE)
  red <- "pca_sub"
}

dims <- seq_len(opt$n_dims)
sub <- RunUMAP(sub, reduction = red, dims = dims,
               reduction.name = "umap_sub", seed.use = 42, verbose = FALSE)
sub <- FindNeighbors(sub, reduction = red, dims = dims,
                     graph.name = c("sub_nn", "sub_snn"), verbose = FALSE)

# ── Optional resolution sweep ─────────────────────────────────────────────────
if (opt$run_sweep) {
  suppressPackageStartupMessages(library(clustree))
  for (r in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)) {
    sub <- FindClusters(sub, graph.name = "sub_snn", resolution = r, random.seed = 42,  verbose = FALSE)
    sub[[glue("subres_{r}")]] <- Idents(sub)
  }
  ct <- clustree(sub@meta.data, prefix = "subres_") + ggtitle("Endothelial subset sweep")
  ggsave(out("09_endo_clustree.pdf"), ct, width = 9, height = 12)
  message(glue("[subset] clustree sweep -> {out('09_endo_clustree.pdf')}"))
}

# ── Cluster at locked resolution ──────────────────────────────────────────────
sub <- FindClusters(sub, graph.name = "sub_snn", resolution = opt$resolution, random.seed = 42, verbose = FALSE)
sub$subcluster <- Idents(sub)
message(glue("[subset] res {opt$resolution} -> {nlevels(sub$subcluster)} subclusters"))

# ── QC overlays: are subclusters beds, or artifacts of batch / cell cycle? ────
p_sub  <- DimPlot(sub, reduction = "umap_sub", group.by = "subcluster",
                  label = TRUE, repel = TRUE) + ggtitle(glue("Subclusters (res {opt$resolution})")) + NoLegend()
p_lin  <- DimPlot(sub, reduction = "umap_sub", group.by = "lineage") + ggtitle(glue("Prior lineage: {paste(lineages, collapse=' + ')}"))
p_time <- DimPlot(sub, reduction = "umap_sub", group.by = "timepoint", cols = tp_cols) +
  ggtitle("Timepoint (trajectory axis / integration QC)")
p_phase <- DimPlot(sub, reduction = "umap_sub", group.by = "Phase") + ggtitle("Cell-cycle phase (artifact check)")
ggsave(out("09_endo_umap_overlays.pdf"), (p_sub | p_lin) / (p_time | p_phase), width = 14, height = 12)

# ── Endothelial bed panel (grouped dotplot; absent genes dropped + reported) ──
bed_panel <- list(
  Progenitor   = c("Etv2", "Tal1", "Lmo2", "Fli1", "Kdr"),
  Pan_EC       = c("Pecam1", "Cdh5", "Tie1", "Cldn5", "Emcn"),
  Arterial     = c("Gja5", "Gja4", "Efnb2", "Dll4", "Hey1", "Sox17"),
  Venous       = c("Nr2f2", "Nrp2", "Aplnr", "Ephb4"),
  Endocardium  = c("Npr3", "Nfatc1", "Gata5", "Klf2", "Bmx"),
  Hemogenic    = c("Runx1", "Gfi1", "Itga2b", "Spi1"),     # EHT -> blood branch
  Lymphatic    = c("Prox1", "Lyve1", "Flt4"),              # expect sparse at E<=10.5
  SHF          = c("Isl1", "Tbx1", "Mef2c", "Fgf8", "Fgf10", "Six2", "Hand2"),
  Pharyngeal_EC = c("Cxcr4","Plxnd1","Gbx2"),
  Integrin_a5b1 = c("Itga5","Itgb1","Fn1"),
  Integrin_aV   = c("Itgav","Itgb3","Itgb5"),
  Integrin_other= c("Itga1","Itga3","Itga4","Itga6","Itga9","Itga8")
)
present   <- rownames(sub)
bed_keep  <- map(bed_panel, ~ .x[.x %in% present])
bed_drop  <- keep(map2(bed_panel, bed_keep, setdiff), ~ length(.x) > 0)
if (length(bed_drop) > 0) {
  message("[panel] bed genes not in object (dropped):")
  iwalk(bed_drop, ~ message(glue("  {.y}: {paste(.x, collapse=', ')}")))
}
bed_keep <- keep(bed_keep, ~ length(.x) > 0)

p_dot <- DotPlot(sub, features = bed_keep, group.by = "subcluster") + RotatedAxis() +
  labs(title = glue("Endothelial bed markers by subcluster (res {opt$resolution})")) +
  theme(axis.text.x = element_text(size = 7), strip.text.x = element_text(size = 7, angle = 90))
ggsave(out("09_endo_bed_dotplot.pdf"), p_dot, width = 18, height = 6)

# ── Unbiased subcluster markers (annotation aid) ──────────────────────────────
future::plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)
Idents(sub) <- "subcluster"
markers <- FindAllMarkers(sub, layer = "data", only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write_csv(markers, out("09_endo_subcluster_markers.csv"))
markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 8, with_ties = FALSE) |>
  ungroup() |> write_csv(out("09_endo_subcluster_top.csv"))

# ── Subcluster composition (timepoint + prior lineage) ────────────────────────
sub@meta.data |> as_tibble() |>
  count(subcluster, timepoint, lineage) |>
  write_csv(out("09_endo_subcluster_summary.csv"))

saveRDS(sub, out("09_endo_subset.rds"))
message(glue("[subset] done -> {out('09_endo_subset.rds')}"))
