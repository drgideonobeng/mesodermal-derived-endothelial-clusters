#!/usr/bin/env Rscript
# =============================================================================
# paa-dev-trajectory — 02: fastMNN integration (temporal merge order)
# -----------------------------------------------------------------------------
# Corrects the per-capture batch the baseline exposed. Runs batchelor::fastMNN
# on log-normalized HVG with merge.order = E80 -> E825 -> E95 -> E105 so each
# step joins developmentally adjacent captures (most shared cell states). The
# corrected embedding is imported into the Seurat object as reduction "mnn".
#
# Re-computes the timepoint mixing matrix on the corrected space and prints a
# before -> after comparison: adjacent-stage overlap should RISE, distant-stage
# overlap should stay near zero. Everything mixing = over-correction.
# =============================================================================
if (!requireNamespace("batchelor", quietly=TRUE))
  stop("batchelor not installed. Run: BiocManager::install('batchelor')")
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(glue)
  library(BiocNeighbors); library(optparse)
  library(batchelor); library(SingleCellExperiment)
})

option_list <- list(
  make_option("--obj", type="character",
    default="~/bioinformatics/scRNAseq/paa-dev-trajectory/results/objects/01_merged_baseline.rds"),
  make_option("--outdir", type="character",
    default="~/bioinformatics/scRNAseq/paa-dev-trajectory/results"),
  make_option("--n_dims", type="integer", default=30L, help="fastMNN corrected dims (match baseline n_pcs)"),
  make_option("--k", type="integer", default=30L),
  make_option("--seed", type="integer", default=42L)
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)

outdir   <- path.expand(opt$outdir)
obj_dir  <- file.path(outdir,"objects"); plot_dir <- file.path(outdir,"plots")
tp_order <- c("E80","E825","E95","E105")

mixing_matrix <- function(emb, tp, levs, k) {
  knn <- BiocNeighbors::findKNN(emb, k=k, BNPARAM=KmknnParam())$index
  m <- matrix(0, length(levs), length(levs), dimnames=list(focal=levs, neighbour=levs))
  for (l in levs) {
    cells <- which(tp==l); if(!length(cells)) next
    nbr <- matrix(tp[knn[cells,,drop=FALSE]], nrow=length(cells))
    for (n in levs) m[l,n] <- mean(rowMeans(nbr==n))
  }
  m
}

message("[02] Loading baseline object...")
obj <- readRDS(path.expand(opt$obj))
obj$timepoint <- factor(obj$timepoint, levels=tp_order)

message("[02] fastMNN (temporal merge order)...")
logmat <- GetAssayData(obj, assay="RNA", layer="data")
hvg    <- VariableFeatures(obj)
batch  <- as.character(obj$timepoint)
set.seed(opt$seed)
mnn <- batchelor::fastMNN(logmat, batch=batch, subset.row=hvg,
                          d=opt$n_dims, merge.order=tp_order)
corrected <- SingleCellExperiment::reducedDim(mnn, "corrected")
if (!is.null(rownames(corrected))) corrected <- corrected[colnames(obj), , drop=FALSE]
colnames(corrected) <- paste0("MNN_", seq_len(ncol(corrected)))
obj[["mnn"]] <- CreateDimReducObject(embeddings=corrected, key="MNN_", assay="RNA")
message(glue("[02] corrected embedding: {ncol(corrected)} dims"))

obj <- RunUMAP(obj, reduction="mnn", dims=1:opt$n_dims,
               reduction.name="umap.mnn", seed.use=opt$seed, verbose=FALSE)

tp_cols <- c(E80="#440154", E825="#31688E", E95="#35B779", E105="#FDE725")
p <- DimPlot(obj, reduction="umap.mnn", group.by="timepoint", cols=tp_cols, pt.size=0.3) +
  ggtitle("fastMNN (temporal merge) - coloured by timepoint")
ggsave(file.path(plot_dir,"02_mnn_umap_timepoint.pdf"), p, width=8, height=7)

tp <- as.character(obj$timepoint)
mix_mnn <- mixing_matrix(corrected, tp, tp_order, opt$k)
message("[02] fastMNN timepoint neighbour-composition:")
print(round(mix_mnn,3))
write.csv(round(mix_mnn,4), file.path(outdir,"02_timepoint_mixing_mnn.csv"))

mix_base <- as.matrix(read.csv(file.path(outdir,"01_timepoint_mixing_baseline.csv"), row.names=1))
dimnames(mix_base) <- list(focal=tp_order, neighbour=tp_order)
adj  <- list(c("E80","E825"), c("E825","E95"), c("E95","E105"))
dist <- list(c("E80","E95"),  c("E80","E105"), c("E825","E105"))
pair_mean <- function(m, prs) mean(sapply(prs, function(p)(m[p[1],p[2]]+m[p[2],p[1]])/2))

message("\n[02] ===== Adjacent vs distant overlap (baseline -> fastMNN) =====")
message(glue("  Adjacent (want INCREASE): {sprintf('%.3f',pair_mean(mix_base,adj))} -> {sprintf('%.3f',pair_mean(mix_mnn,adj))}"))
message(glue("  Distant  (want ~0):       {sprintf('%.3f',pair_mean(mix_base,dist))} -> {sprintf('%.3f',pair_mean(mix_mnn,dist))}"))

mix_df <- as.data.frame(as.table(mix_mnn)); names(mix_df) <- c("focal","neighbour","frac")
mix_df$focal <- factor(mix_df$focal, levels=tp_order)
mix_df$neighbour <- factor(mix_df$neighbour, levels=tp_order)
pm <- ggplot(mix_df, aes(neighbour, focal, fill=frac)) +
  geom_tile() + geom_text(aes(label=sprintf("%.2f",frac)), size=3.5) +
  scale_fill_viridis_c(limits=c(0,1), name="kNN frac") +
  labs(title=glue("Timepoint neighbour composition (k={opt$k}) - fastMNN"),
       subtitle="Adjacent stages should now mix; distant stages should stay apart",
       x="neighbour timepoint", y="focal timepoint") +
  theme_minimal(base_size=11) + coord_equal()
ggsave(file.path(plot_dir,"02_timepoint_mixing_mnn.pdf"), pm, width=6.5, height=5.5)

saveRDS(obj, file.path(obj_dir,"02_merged_mnn.rds"))
message(glue("\n[02] Done -> {file.path(obj_dir,'02_merged_mnn.rds')}"))
