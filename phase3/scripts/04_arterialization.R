#!/usr/bin/env Rscript
# 04_arterialization.R
# Objective 2: the late-window plexus -> artery remodeling gene program.
# Continuous ordering of the late developmental window by arterial identity,
# rooted at the LEAST-arterialized cells (an output of this script, not an
# assumption), with monotonically-changing genes extracted as the remodeling
# program candidate set.
#
# ── WHY CONTINUOUS ORDERING, NOT CLUSTER-CHAINING (read before running) ───────
# Slingshot on the full uncorrected embedding (03_trajectory.R) gave a
# biologically implausible topology (arterial-from-venous) because (a) stage
# gaps fragment beds across UMAP islands and (b) EC identity is a continuum --
# clustering manufactures discreteness that isn't structurally there. Rather
# than coax a graph-based trajectory into the late window, this script orders
# cells directly by a continuous, integrin-free ARTERIAL IDENTITY SCORE. This
# is content-based ordering (what the cell expresses), not geometry-based
# (where the cell sits in a kNN graph) -- it sidesteps the fragmentation
# problem entirely. PCA/UMAP below are for QC VISUALIZATION ONLY; they are
# NOT the ordering axis.
#
# ── VALIDATING THE SCORE AXIS (new) ────────────────────────────────────────────
# Choosing a curated score over the window's own PCA is a real assumption:
# nothing guarantees arterialization is the dominant axis of variance even
# after narrowing to two adjacent timepoints (residual stage/technical
# difference and cell cycle are still plausible competitors). Section 3a below
# checks this directly -- correlating Arterial_score against the late-window's
# own top PCs. Strong correlation is independent evidence the score tracks
# real dominant structure, not just an assumed one; weak correlation flags
# something else may be driving this window's variance instead (see
# 04_score_vs_pca_diagnostic.*).
#
# ── WHY THIS WINDOW, WHY NO RE-INTEGRATION ─────────────────────────────────────
# Late window = native uncorrected clusters {0, 5} (E9.5/E10.5; verified in
# 03a_traj_native_markers.R). Restricting to two adjacent, already-continuous
# timepoints is what makes "continuous" a defensible description in the first
# place. We deliberately do NOT fastMNN-integrate here: with only two
# timepoints left, correction risks flattening the very within-window
# variation this script exists to characterize (the same over-correction
# concern documented in 01_subset_endothelium.R).
#
# ── WHY DROP Hb GENES FROM HVGs ────────────────────────────────────────────────
# Native cluster 0 carries severe ambient hemoglobin contamination (Hba/Hbb
# 68-86% of cells positive; diagnosed during 03a verification). Left in, Hb
# genes can dominate variable-feature selection and drive the embedding/score
# on an artifact rather than biology. They are excluded from HVGs (and from
# the arterial score candidate pool by construction); raw counts are
# untouched, so per-cluster Hb diagnostics below still work.
#
# ── WHY ROOT BY SCORE, NOT BY ASSUMED IDENTITY ─────────────────────────────────
# Venous origin of the low-arterial population was already tested and came
# back null (Hypothesis_testing_SHF_v_nSHF.R: venous score LOWEST, not
# highest, in SHF-low arterial cells). So the root here is defined purely as
# "--root_quantile of cells with the lowest arterial score," and its identity
# is then CHARACTERIZED post hoc via FindMarkers + composition tables --
# an output, not a label applied in advance.
#
# ── ANTI-DOUBLE-DIPPING ────────────────────────────────────────────────────────
# Genes used to BUILD the arterial score (--arterial_markers) will trivially
# correlate with it -- that's a sanity check the score worked, not a
# discovery. They stay in the full correlation table (flagged via
# in_score_panel) but are EXCLUDED from the declared "remodeling program"
# gene list, mirroring the same anti-circularity logic in 02_shf_modules.R.
#
# ── CAVEATS ─────────────────────────────────────────────────────────────────────
# • Late window is still drawn from two discrete native clusters. Check the
#   arterial-score histogram (04_arterial_score_distribution.pdf) for a
#   continuum vs bimodality before trusting "continuous" -- if it's bimodal,
#   the window may not be a single remodeling process.
# • Small n: late-window cell counts are modest, so correlation power is
#   limited. Treat p_adj as a screen, not a verdict -- inspect effect size
#   (rho) and pct_expressing together.
# • Cell-level Hb-high filtering (beyond HVG exclusion) is NOT applied here;
#   --hb_pct_filter is provided but off by default. Turn on only if ambient
#   Hb is still visibly distorting the late-window UMAP after HVG exclusion.
#
# Input : results/phase3/<lineage>/uncorrected/01_endo_subset.rds
# Output: 04_late_window_subset.rds | 04_arterial_score_distribution.pdf |
#         04_umap_overlays.pdf | 04_root_markers.csv |
#         04_gene_arterialization_correlations.csv |
#         04_score_vs_pca_diagnostic.csv | 04_score_vs_pca_diagnostic.pdf |
#         04_top_monotonic_genes.pdf | 04_remodeling_summary.txt

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage",         type = "character", default = "endothelium",
              help = "phase-3 namespace under results/phase3/<lineage>/"),
  make_option("--input",           type = "character", default = NULL,
              help = "uncorrected subset; default results/phase3/<lineage>/uncorrected/01_endo_subset.rds"),
  make_option("--outdir",          type = "character", default = NULL,
              help = "default results/phase3/<lineage>/uncorrected/late_window"),
  make_option("--late_clusters",   type = "character", default = "0,5",
              help = "comma-separated native subcluster IDs defining the late window (E9.5/E10.5)"),
  make_option("--hb_genes",        type = "character",
              default = "Hba-a1,Hba-a2,Hba-x,Hbb-bs,Hbb-bt,Hbb-y,Hbb-bh1",
              help = "mouse hemoglobin gene symbols dropped from HVG selection"),
  make_option("--hb_pct_filter",   type = "double", default = NULL,
              help = "OFF by default; if set (0-1), drop cells with >this fraction of UMIs from --hb_genes"),
  make_option("--arterial_markers", type = "character",
              default = "Gja5,Tmem100,Sox17,Dll4,Hey1,Efnb2,Cxcr4,Epas1",
              help = "integrin-free arterial identity panel (module score + ordering axis)"),
  make_option("--n_hvgs",          type = "integer", default = 1500L,
              help = "HVGs re-selected within the late window (smaller n than 01: fewer cells)"),
  make_option("--n_dims",          type = "integer", default = 15L,
              help = "PCA/UMAP dims for QC visualization only -- not used for ordering"),
  make_option("--root_quantile",   type = "double", default = 0.2,
              help = "bottom quantile of arterial score defining the 'root' group for FindMarkers"),
  make_option("--min_pct",         type = "double", default = 0.1,
              help = "FindMarkers min.pct for root-state characterization"),
  make_option("--logfc_threshold", type = "double", default = 0.25,
              help = "FindMarkers logFC threshold for root-state characterization"),
  make_option("--cor_method",      type = "character", default = "spearman"),
  make_option("--min_rho",         type = "double", default = 0.3,
              help = "|rho| threshold for the declared monotonic remodeling-program gene list"),
  make_option("--alpha",           type = "double", default = 0.05,
              help = "BH-adjusted p threshold for the declared monotonic gene list")
)))

set.seed(42)

# ── Resolve paths from --lineage ────────────────────────────────────────────────
base <- file.path("results/phase3", opt$lineage)
if (is.null(opt$input))  opt$input  <- file.path(base, "uncorrected", "01_endo_subset.rds")
if (is.null(opt$outdir)) opt$outdir <- file.path(base, "uncorrected", "late_window")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

late_clusters <- str_split(opt$late_clusters, ",")[[1]] |> str_trim()
hb_genes      <- str_split(opt$hb_genes, ",")[[1]] |> str_trim()
arterial_genes_req <- str_split(opt$arterial_markers, ",")[[1]] |> str_trim()

# ══════════════════════════════════════════════════════════════════════════════
# 1. SUBSET TO THE LATE WINDOW
# ══════════════════════════════════════════════════════════════════════════════
obj <- readRDS(opt$input)
stopifnot("subcluster" %in% colnames(obj@meta.data))

miss_clu <- setdiff(late_clusters, levels(obj$subcluster))
if (length(miss_clu) > 0)
  stop(glue("[window] requested late cluster(s) not present: {paste(miss_clu, collapse=', ')} -- ",
            "have: {paste(levels(obj$subcluster), collapse=', ')}"))

late <- subset(obj, subset = subcluster %in% late_clusters)
late$prior_subcluster <- droplevels(late$subcluster)
message(glue("[window] {ncol(late)} cells in late window (clusters {paste(late_clusters, collapse=' + ')})"))
print(table(subcluster = late$prior_subcluster, timepoint = late$timepoint))

# ── Hb contamination diagnostic (by native cluster, pre-exclusion) ────────────
hb_present <- intersect(hb_genes, rownames(late))
if (length(hb_present) > 0) {
  hb_pct <- FetchData(late, c(hb_present, "prior_subcluster")) |>
    as_tibble() |>
    pivot_longer(-prior_subcluster, names_to = "gene", values_to = "expr") |>
    group_by(prior_subcluster, gene) |>
    summarise(pct_pos = mean(expr > 0) * 100, .groups = "drop")
  message("[Hb] %cells+ by native cluster (ambient contamination check):")
  print(hb_pct |> pivot_wider(names_from = gene, values_from = pct_pos))
} else {
  message("[Hb] none of --hb_genes found in object -- skipping diagnostic")
}

# ── Optional cell-level Hb filter (off by default; see header caveat) ─────────
if (!is.null(opt$hb_pct_filter) && length(hb_present) > 0) {
  hb_frac <- Matrix::colSums(GetAssayData(late, layer = "counts")[hb_present, , drop = FALSE]) /
             late$nCount_RNA
  keep <- hb_frac <= opt$hb_pct_filter
  message(glue("[Hb] --hb_pct_filter {opt$hb_pct_filter}: dropping {sum(!keep)} / {ncol(late)} cells"))
  late <- subset(late, cells = colnames(late)[keep])
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. RE-NORMALIZE + RE-SELECT HVGs WITHIN THE WINDOW (Hb genes excluded)
# ══════════════════════════════════════════════════════════════════════════════
DefaultAssay(late) <- "RNA"
late[["RNA"]] <- JoinLayers(late[["RNA"]])
late <- late |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = opt$n_hvgs, verbose = FALSE)

vf <- setdiff(VariableFeatures(late), hb_present)
VariableFeatures(late) <- vf
message(glue("[HVG] {length(vf)} variable features after excluding {length(hb_present)} Hb gene(s)"))

# ── Fresh PCA + UMAP -- QC VISUALIZATION ONLY, not the ordering axis ──────────
late <- late |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = opt$n_dims, reduction.name = "pca_late", verbose = FALSE) |>
  RunUMAP(reduction = "pca_late", dims = seq_len(opt$n_dims),
          reduction.name = "umap_late", seed.use = 42, verbose = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 3. ARTERIAL IDENTITY SCORE -- the ordering axis
# ══════════════════════════════════════════════════════════════════════════════
arterial_genes <- arterial_genes_req[arterial_genes_req %in% rownames(late)]
if (length(setdiff(arterial_genes_req, arterial_genes)) > 0)
  message(glue("[score] arterial panel gene(s) not in object (dropped): ",
               "{paste(setdiff(arterial_genes_req, arterial_genes), collapse=', ')}"))

late <- AddModuleScore(late, features = list(Arterial = arterial_genes),
                       name = "Art", ctrl = 50, seed = 42)
late$Arterial_score <- late$Art1; late$Art1 <- NULL
late$arterial_rank  <- rank(late$Arterial_score) / ncol(late)   # 0 (least) -> 1 (most arterialized)

message(glue("[score] Arterial score range: [{round(min(late$Arterial_score),3)}, ",
             "{round(max(late$Arterial_score),3)}]; panel: {paste(arterial_genes, collapse=', ')}"))

# ── Root group: bottom --root_quantile of arterial score ──────────────────────
cut <- quantile(late$Arterial_score, opt$root_quantile)
late$window_state <- ifelse(late$Arterial_score <= cut, "Low_arterial_root", "Rest")
message(glue("[root] bottom {opt$root_quantile*100}% (score <= {round(cut,3)}): ",
             "{sum(late$window_state=='Low_arterial_root')} cells"))
print(table(root = late$window_state, timepoint = late$timepoint,
            native_cluster = late$prior_subcluster))

# ══════════════════════════════════════════════════════════════════════════════
# 3a. VALIDATE: does Arterial_score track this window's OWN dominant unsupervised
#     structure (pca_late), or could something else (residual stage/technical
#     difference, cell cycle) be the bigger signal here? Does NOT change the
#     ordering axis -- this is a sanity check on it, not a replacement for it.
# ══════════════════════════════════════════════════════════════════════════════
pca_emb    <- Embeddings(late, "pca_late")
n_pc_check <- min(5L, ncol(pca_emb))

pc_cor <- map_df(seq_len(n_pc_check), function(i) {
  ct <- suppressWarnings(cor.test(pca_emb[, i], late$Arterial_score, method = opt$cor_method))
  tibble(PC = glue("PC{i}"), rho = unname(ct$estimate), p = ct$p.value)
}) |> mutate(p_adj = p.adjust(p, "BH"))
write_csv(pc_cor, out("04_score_vs_pca_diagnostic.csv"))

best_pc <- pc_cor |> slice_max(abs(rho), n = 1)
message(glue("[validate] Arterial score vs late-window PCs (top {n_pc_check}):"))
print(pc_cor)
if (abs(best_pc$rho) >= 0.5) {
  message(glue("[validate] {best_pc$PC} correlates with Arterial_score (rho={round(best_pc$rho,2)}) -- ",
               "arterialization looks like a real component of this window's dominant structure."))
} else {
  message(glue("[validate] WEAK correlation with every top PC (best: {best_pc$PC}, rho={round(best_pc$rho,2)}) -- ",
               "the window's biggest unsupervised signal may be something other than arterialization ",
               "(residual stage/technical difference, cell cycle). Inspect 04_score_vs_pca_diagnostic.pdf ",
               "and consider checking Phase / nCount_RNA against the top PC before trusting the score-based ordering."))
}

pc_plot_df <- map_df(seq_len(n_pc_check), function(i) {
  tibble(PC = glue("PC{i}"), pc_value = pca_emb[, i], Arterial_score = late$Arterial_score)
}) |> mutate(PC = factor(PC, levels = glue("PC{seq_len(n_pc_check)}")))

p_pcval <- ggplot(pc_plot_df, aes(pc_value, Arterial_score)) +
  geom_point(size = .4, alpha = .4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6, colour = "firebrick") +
  facet_wrap(~ PC, scales = "free_x", nrow = 1) +
  labs(title = "Arterial score vs late-window PCA (validation, not the ordering axis)",
       subtitle = "does the curated score track this window's own dominant unsupervised structure?",
       x = "PC value", y = "Arterial score") +
  theme_minimal(base_size = 10)
ggsave(out("04_score_vs_pca_diagnostic.pdf"), p_pcval, width = 3 * n_pc_check, height = 4)

# ══════════════════════════════════════════════════════════════════════════════
# 4. CHARACTERIZE THE ROOT STATE (output, not assumption)
# ══════════════════════════════════════════════════════════════════════════════
Idents(late) <- "window_state"
root_markers <- FindMarkers(late, ident.1 = "Low_arterial_root", ident.2 = "Rest",
                            logfc.threshold = opt$logfc_threshold, min.pct = opt$min_pct,
                            verbose = FALSE) |>
  rownames_to_column("gene") |> arrange(desc(avg_log2FC)) |> as_tibble()
write_csv(root_markers, out("04_root_markers.csv"))
message(glue("[root] top markers (positive logFC, root vs rest):"))
print(head(root_markers |> filter(avg_log2FC > 0), 15))

# ══════════════════════════════════════════════════════════════════════════════
# 5. GENE x ARTERIAL-SCORE CORRELATION -- the remodeling program candidates
#    Candidate universe = window-specific HVGs (genes variable WITHIN this
#    window -- the natural candidate pool for a within-window program).
# ══════════════════════════════════════════════════════════════════════════════
candidate_genes <- union(VariableFeatures(late), arterial_genes)
expr_mat <- as.matrix(GetAssayData(late, layer = "data")[candidate_genes, , drop = FALSE])
score_vec <- late$Arterial_score

cor_res <- map_df(candidate_genes, function(g) {
  x  <- expr_mat[g, ]
  ct <- suppressWarnings(cor.test(x, score_vec, method = opt$cor_method))
  tibble(gene = g, rho = unname(ct$estimate), p = ct$p.value,
         pct_expressing = mean(x > 0))
}) |>
  mutate(p_adj = p.adjust(p, "BH"),
         in_score_panel = gene %in% arterial_genes) |>
  arrange(desc(abs(rho)))
write_csv(cor_res, out("04_gene_arterialization_correlations.csv"))

# ── Declared remodeling-program gene list: significant, |rho| >= threshold,
#    EXCLUDING score-panel genes (anti-double-dipping) ────────────────────────
program <- cor_res |>
  filter(!in_score_panel, p_adj < opt$alpha, abs(rho) >= opt$min_rho) |>
  arrange(desc(abs(rho)))
message(glue("[program] {nrow(program)} gene(s) pass |rho|>={opt$min_rho}, p_adj<{opt$alpha} ",
             "(score-panel genes excluded by design)"))
print(head(program, 20))

# ══════════════════════════════════════════════════════════════════════════════
# 6. PLOTS
# ══════════════════════════════════════════════════════════════════════════════
# -- Arterial score distribution (continuum vs bimodal QC check) --------------
p_dist <- ggplot(late@meta.data, aes(Arterial_score, fill = timepoint)) +
  geom_histogram(bins = 40, alpha = .8, position = "stack") +
  geom_vline(xintercept = cut, linetype = "dashed", colour = "grey20") +
  labs(title = "Arterial score distribution -- late window",
       subtitle = glue("dashed line = root cutoff (bottom {opt$root_quantile*100}%); ",
                        "inspect for continuum vs bimodality"),
       x = "Arterial score", y = "Cells") +
  theme_minimal(base_size = 11)
ggsave(out("04_arterial_score_distribution.pdf"), p_dist, width = 8, height = 5)

# -- UMAP overlays: score / timepoint / root status / native cluster ----------
p1 <- FeaturePlot(late, "Arterial_score", reduction = "umap_late") +
  ggtitle("Arterial score")
p2 <- DimPlot(late, reduction = "umap_late", group.by = "timepoint") +
  ggtitle("Timepoint")
p3 <- DimPlot(late, reduction = "umap_late", group.by = "window_state") +
  ggtitle("Root status")
p4 <- DimPlot(late, reduction = "umap_late", group.by = "prior_subcluster") +
  ggtitle("Native cluster (prior)")
ggsave(out("04_umap_overlays.pdf"), (p1 | p2) / (p3 | p4), width = 12, height = 10)

# -- Top monotonic genes vs arterial rank (declared program only) -------------
if (nrow(program) > 0) {
  top_genes <- head(program$gene, min(20, nrow(program)))
  plot_df <- FetchData(late, c(top_genes, "arterial_rank")) |>
    as_tibble() |>
    pivot_longer(-arterial_rank, names_to = "gene", values_to = "expr") |>
    mutate(gene = factor(gene, levels = top_genes))
  p_mono <- ggplot(plot_df, aes(arterial_rank, expr)) +
    geom_point(size = .3, alpha = .4) +
    geom_smooth(method = "loess", se = FALSE, linewidth = .6, colour = "firebrick") +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(title = "Top monotonic genes along the arterialization axis",
         subtitle = "score-panel genes excluded; x = arterial rank (0=least, 1=most arterialized)",
         x = "Arterial rank", y = "Expression") +
    theme_minimal(base_size = 9)
  ggsave(out("04_top_monotonic_genes.pdf"), p_mono, width = 12, height = 10)
} else {
  message("[plot] no genes passed the program threshold -- skipping 04_top_monotonic_genes.pdf")
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. SUMMARY + SAVE
# ══════════════════════════════════════════════════════════════════════════════
writeLines(c(
  glue("Late window: {ncol(late)} cells, native clusters {paste(late_clusters, collapse=' + ')}"),
  glue("Arterial panel ({length(arterial_genes)} genes): {paste(arterial_genes, collapse=', ')}"),
  glue("Root group (bottom {opt$root_quantile*100}% arterial score): ",
       "{sum(late$window_state=='Low_arterial_root')} cells"),
  glue("Root state top positive markers: ",
       "{paste(head(root_markers$gene[root_markers$avg_log2FC>0], 10), collapse=', ')}"),
  glue("Score-vs-PCA validation: best correlate {best_pc$PC} (rho={round(best_pc$rho,2)}, p_adj={signif(best_pc$p_adj,3)})"),
  glue("Remodeling program (|rho|>={opt$min_rho}, p_adj<{opt$alpha}, panel genes excluded): ",
       "{nrow(program)} gene(s)"),
  glue("Top program genes: {paste(head(program$gene, 15), collapse=', ')}")
), out("04_remodeling_summary.txt"))

saveRDS(late, out("04_late_window_subset.rds"))
message(glue("[done] -> {out('04_late_window_subset.rds')}"))