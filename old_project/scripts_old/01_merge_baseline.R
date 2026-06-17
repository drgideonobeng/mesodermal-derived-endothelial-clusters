#!/usr/bin/env Rscript
# =============================================================================
# paa-dev-trajectory — 01: merge + no-integration baseline
# -----------------------------------------------------------------------------
# Loads the four GSE158941 WT timecourse samples (post-QC objects from
# paa-development Phase 1), merges them, and builds a NO-INTEGRATION baseline
# embedding: LogNormalize -> HVG -> PCA -> UMAP.
#
# Also computes a timepoint k-NN "mixing matrix": for cells of each timepoint,
# the mean fraction of their nearest neighbours from each timepoint (in PCA
# space). Expectation for a timecourse: adjacent stages (E80/E825) mix heavily,
# distant stages (E80/E105) stay apart. That graded structure is the yardstick
# for judging fastMNN in script 02 — and tells us now whether no-integration
# already suffices. A near-zero E80<->E825 overlap would flag a technical split
# that fastMNN should fix.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
  library(glue);   library(BiocNeighbors); library(optparse)
})

option_list <- list(
  make_option("--src_dir", type="character",
    default="~/bioinformatics/scRNAseq/paa-development/results/phase1",
    help="Phase 1 dir with <sample>/objects/04_seurat_filtered.rds"),
  make_option("--outdir",  type="character",
    default="~/bioinformatics/scRNAseq/paa-dev-trajectory/results"),
  make_option("--n_hvg", type="integer", default=2000L),
  make_option("--n_pcs", type="integer", default=30L),
  make_option("--k",     type="integer", default=30L, help="k for mixing diagnostic"),
  make_option("--seed",  type="integer", default=42L)
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)

src      <- path.expand(opt$src_dir)
outdir   <- path.expand(opt$outdir)
obj_dir  <- file.path(outdir, "objects"); plot_dir <- file.path(outdir, "plots")
dir.create(obj_dir,  recursive=TRUE, showWarnings=FALSE)
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

# Four GSE158941 WT samples in developmental order
tp_order <- c("E80","E825","E95","E105")
samples <- tibble::tribble(
  ~sample_id,       ~timepoint,
  "mesp1_e80",      "E80",
  "mesp1_e825",     "E825",
  "mesp1_e95_rep1", "E95",
  "mesp1_e105",     "E105"
)

message("[01] Loading four GSE158941 objects...")
objs <- lapply(seq_len(nrow(samples)), function(i) {
  f <- file.path(src, samples$sample_id[i], "objects", "04_seurat_filtered.rds")
  stopifnot(file.exists(f))
  o <- readRDS(f); DefaultAssay(o) <- "RNA"
  o$sample_id <- samples$sample_id[i]; o$timepoint <- samples$timepoint[i]
  message(glue("[01]   {samples$sample_id[i]}: {ncol(o)} cells")); o
})

merged <- merge(objs[[1]], y=objs[-1], add.cell.ids=samples$sample_id, project="paa_traj")
DefaultAssay(merged) <- "RNA"
if (packageVersion("Seurat") >= "5.0.0") merged <- JoinLayers(merged)
merged$timepoint <- factor(merged$timepoint, levels=tp_order)

message(glue("[01] Merged: {ncol(merged)} cells x {nrow(merged)} genes"))
merged@meta.data |>
  group_by(timepoint) |>
  summarise(n=n(),
            median_nFeature=median(nFeature_RNA),
            median_nCount=median(nCount_RNA), .groups="drop") |>
  as.data.frame() |> print()

message("[01] LogNormalize -> HVG -> scale -> PCA -> UMAP (no integration)...")
merged <- merged |>
  NormalizeData(verbose=FALSE) |>
  FindVariableFeatures(nfeatures=opt$n_hvg, verbose=FALSE) |>
  ScaleData(verbose=FALSE) |>
  RunPCA(npcs=opt$n_pcs, verbose=FALSE)
merged <- RunUMAP(merged, dims=1:opt$n_pcs, seed.use=opt$seed, verbose=FALSE)

tp_cols <- c(E80="#440154", E825="#31688E", E95="#35B779", E105="#FDE725")
p_umap <- DimPlot(merged, group.by="timepoint", cols=tp_cols, pt.size=0.3) +
  ggtitle("No-integration baseline - coloured by timepoint")
ggsave(file.path(plot_dir, "01_baseline_umap_timepoint.pdf"), p_umap, width=8, height=7)

# ---- Timepoint k-NN mixing matrix (PCA space) -------------------------------
message(glue("[01] Mixing matrix (k={opt$k}, PCA space)..."))
emb <- Embeddings(merged, "pca")[, 1:opt$n_pcs]
knn <- BiocNeighbors::findKNN(emb, k=opt$k, BNPARAM=KmknnParam())$index
tp  <- as.character(merged$timepoint); levs <- tp_order
mix <- matrix(0, length(levs), length(levs), dimnames=list(focal=levs, neighbour=levs))
for (l in levs) {
  cells <- which(tp == l); if (!length(cells)) next
  nbr <- matrix(tp[knn[cells, , drop=FALSE]], nrow=length(cells))
  for (m in levs) mix[l, m] <- mean(rowMeans(nbr == m))
}
message("[01] Timepoint neighbour-composition (rows sum to 1):")
print(round(mix, 3))
write.csv(round(mix, 4), file.path(outdir, "01_timepoint_mixing_baseline.csv"))

mix_df <- as.data.frame(as.table(mix)); names(mix_df) <- c("focal","neighbour","frac")
mix_df$focal     <- factor(mix_df$focal, levels=tp_order)
mix_df$neighbour <- factor(mix_df$neighbour, levels=tp_order)
p_mix <- ggplot(mix_df, aes(neighbour, focal, fill=frac)) +
  geom_tile() + geom_text(aes(label=sprintf("%.2f", frac)), size=3.5) +
  scale_fill_viridis_c(limits=c(0,1), name="kNN frac") +
  labs(title=glue("Timepoint neighbour composition (k={opt$k}) - baseline"),
       subtitle="Adjacent stages should mix; distant stages should not",
       x="neighbour timepoint", y="focal timepoint") +
  theme_minimal(base_size=11) + coord_equal()
ggsave(file.path(plot_dir, "01_timepoint_mixing_baseline.pdf"), p_mix, width=6.5, height=5.5)

saveRDS(merged, file.path(obj_dir, "01_merged_baseline.rds"))
message(glue("[01] Done -> {file.path(obj_dir,'01_merged_baseline.rds')}"))
