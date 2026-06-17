#!/usr/bin/env Rscript
# ============================================================================
# 04_cluster_markers.R   (v3)
# ----------------------------------------------------------------------------
# STEP 2 of whole-manifold annotation: WITHIN-DATA evidence per cluster.
#   (A) Cell-cycle scoring -- CellCycleScoring (Tirosh S/G2M gene sets,
#       mouse-converted). Cell-cycle is a STATE not an identity, so it is
#       handled separately from the lineage panel, which contains only
#       identity markers.
#   (B) Unbiased DE   -- FindAllMarkers (what actually distinguishes each
#       cluster, no prior assumptions)
#   (C) Canonical map -- broad curated panel: grouped dotplot + per-lineage
#       module-score heatmap. v3 changes: Proliferation removed; Venous_EC
#       and Lymphatic_EC are now separate; Notochord_axial carries both T
#       and Tbxt (whichever annotation the object uses is retained);
#       dropped genes are reported explicitly.
#
# Inputs : results/objects/03_clustered.rds
# Outputs: results/04_markers_all.csv / 04_markers_top.csv
#          results/04_lineage_scores_by_cluster.csv
#          results/04_cellcycle_by_cluster.csv
#          results/plots/04_umap_phase.pdf
#          results/plots/04_cluster_phase_composition.pdf
#          results/plots/04_dotplot_datadriven.pdf
#          results/plots/04_dotplot_canonical.pdf
#          results/plots/04_lineage_score_heatmap.pdf
#
# Usage:
#   conda run -n paa-dev Rscript scripts/04_cluster_markers.R \
#       --input results/objects/03_clustered.rds --resolution 0.4
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(glue)
  library(optparse)
  library(Matrix)
})

option_list <- list(
  make_option("--input",      type = "character",
              default = "results/objects/03_clustered.rds"),
  make_option("--outdir",     type = "character", default = "results"),
  make_option("--resolution", type = "double",    default = 0.4),
  make_option("--assay",      type = "character", default = "RNA"),
  make_option("--top_n",      type = "integer",   default = 5L),
  make_option("--min_pct",    type = "double",    default = 0.25),
  make_option("--logfc",      type = "double",    default = 0.25)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir_plot <- file.path(opt$outdir, "plots")
dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# LINEAGE PANEL  (identity markers only -- no cell-cycle genes)
#
# Changes from v2:
#   - Proliferation REMOVED: cell cycle is a state, not an identity; handled
#     below via CellCycleScoring on the Tirosh gene sets.
#   - Venous_lymph_EC SPLIT into Venous_EC and Lymphatic_EC: they are
#     distinct cell types (lymphatic buds from venous ~E9.5-10.5) and
#     separating them lets us confirm the expected lymphatic null at these
#     stages, rather than blending signals.
#   - Notochord_axial: both "T" and "Tbxt" listed (same gene, annotation
#     varies by reference genome); the script retains whichever is present.
#   - Caveat: Nr2f2 (Venous_EC), Prox1 and Lyve1 (Lymphatic_EC) are not
#     endothelial-exclusive -- these calls are only interpretable within
#     Pecam1/Cdh5-positive clusters. Cross-check against the pan-EC panel
#     and data-driven markers before labelling a cluster venous or lymphatic.
# ----------------------------------------------------------------------------
panel <- list(
  Cardiomyocyte    = c("Tnnt2","Myh6","Actc1","Nkx2-5","Tnni3","Myl7"),
  First_heart_fld  = c("Tbx5","Hand1","Nppa"),
  Second_heart_fld = c("Isl1","Fgf8","Fgf10","Six2","Mef2c"),
  Pharyngeal_meso  = c("Tbx1","Pitx2","Tcf21","Lhx2","Msc"),
  Skeletal_muscle  = c("Myf5","Myod1","Myog","Tnnt3"),
  Paraxial_somite  = c("Meox1","Tcf15","Pax3","Foxc1","Foxc2"),
  Endothelium_pan  = c("Pecam1","Cdh5","Kdr","Tie1","Cd34"),
  Angioblast_prog  = c("Etv2","Tal1","Lmo2","Fli1"),
  Arterial_EC      = c("Efnb2","Dll4","Gja5","Hey1"),
  Venous_EC        = c("Nr2f2","Ephb4","Emcn","Aplnr","Nrp2"),
  Lymphatic_EC     = c("Prox1","Lyve1","Pdpn","Flt4"),
  Endocardium      = c("Npr3","Nfatc1"),
  Erythroid_blood  = c("Hba-x","Hbb-bh1","Hbb-y","Hba-a1","Gata1","Klf1","Runx1","Itga2b"),
  LPM_limb         = c("Hand2","Foxf1","Prrx1","Tbx4"),
  Epicardium       = c("Wt1","Tbx18","Upk3b"),
  Smooth_mural     = c("Acta2","Tagln","Myh11","Pdgfrb"),
  Mesenchyme_EMT   = c("Twist1","Snai1","Postn","Col1a1"),
  Notochord_axial  = c("T","Tbxt","Noto","Shh"),
  Neural_ectoderm  = c("Sox2","Pax6","Tubb3","Sox1"),
  Endoderm         = c("Foxa2","Sox17","Epcam","Cldn6","Afp")
)

# ---- Load ------------------------------------------------------------------
message(glue("[load] {opt$input}"))
obj <- readRDS(opt$input)
DefaultAssay(obj) <- opt$assay
message(glue("[assay] {DefaultAssay(obj)} | class = {class(obj[[opt$assay]])[1]}"))
message(glue("[layers] before join: {paste(Layers(obj[[opt$assay]]), collapse=', ')}"))
obj[[opt$assay]] <- tryCatch(
  JoinLayers(obj[[opt$assay]]),
  error = function(e) { message(glue("[join] {conditionMessage(e)}")); obj[[opt$assay]] }
)
message(glue("[layers] after  join: {paste(Layers(obj[[opt$assay]]), collapse=', ')}"))

if (!"data" %in% Layers(obj[[opt$assay]])) {
  message("[norm] no 'data' layer -> NormalizeData")
  obj <- NormalizeData(obj, verbose = FALSE)
}
dat <- GetAssayData(obj, assay = opt$assay, layer = "data")
message(glue("[check] data {nrow(dat)}x{ncol(dat)} | nonzero={nnzero(dat)} | max={round(max(dat),2)}"))
if (nnzero(dat) == 0) stop("[fatal] RNA 'data' layer is all zeros.")
if (max(dat) > 50 && isTRUE(all(dat@x == floor(dat@x)))) {
  message("[norm] looks like raw counts -> re-normalizing")
  obj <- NormalizeData(obj, verbose = FALSE)
}

# ---- Identity --------------------------------------------------------------
res_col <- glue("clust_res{opt$resolution}")
if (!res_col %in% colnames(obj@meta.data))
  stop(glue("'{res_col}' not found. Available: ",
            "{paste(grep('clust_res', colnames(obj@meta.data), value=TRUE), collapse=', ')}"))
Idents(obj) <- res_col
n_clust <- nlevels(Idents(obj))
message(glue("[identity] {res_col} ({n_clust} clusters)"))

present <- rownames(obj)

# ---- Panel filtering with explicit dropped-gene report ---------------------
panel_raw <- panel
panel     <- map(panel, ~ .x[.x %in% present])
dropped   <- map2(panel_raw, panel, ~ setdiff(.x, .y))
dropped   <- dropped[map_lgl(dropped, ~ length(.x) > 0)]
if (length(dropped) > 0) {
  message("[panel] genes NOT found in object (dropped):")
  iwalk(dropped, ~ message(glue("  {.y}: {paste(.x, collapse=', ')}")))
} else {
  message("[panel] all panel genes found in object")
}
panel <- panel[map_lgl(panel, ~ length(.x) > 0)]
message(glue("[panel] {length(unlist(panel))} genes / {length(panel)} lineages retained"))

# ---- (A) Cell-cycle scoring ------------------------------------------------
# CellCycleScoring uses Tirosh et al. human S-phase and G2M gene sets.
# Convert to mouse symbols via str_to_title (e.g. MKI67 -> Mki67); Seurat
# silently drops any that aren't in the object.
message("[CC] scoring cell cycle phases ...")
s_genes   <- str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
message(glue("[CC] S genes in object  : {sum(s_genes %in% present)}/{length(s_genes)}"))
message(glue("[CC] G2M genes in object: {sum(g2m_genes %in% present)}/{length(g2m_genes)}"))

obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                        set.ident = FALSE)
# Summarise phase distribution per cluster
cc_comp <- obj@meta.data |>
  as_tibble() |>
  count(cluster = .data[[res_col]], Phase) |>
  group_by(cluster) |>
  mutate(frac = n / sum(n)) |>
  ungroup()
write_csv(cc_comp, file.path(opt$outdir, "04_cellcycle_by_cluster.csv"))

phase_cols <- c(G1 = "#4575b4", S = "#d73027", G2M = "#fdae61")

p_phase_umap <- DimPlot(obj, reduction = "umap.mnn", group.by = "Phase",
                        cols = phase_cols) +
  labs(title = "Cell cycle phase (Tirosh gene sets, mouse-converted)")

p_phase_comp <- ggplot(cc_comp, aes(cluster, frac, fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = phase_cols) +
  labs(x = glue("cluster ({res_col})"), y = "fraction of cluster",
       title = "Cell cycle phase composition by cluster",
       subtitle = "High S/G2M in a cluster -> proliferative progenitor state") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(dir_plot, "04_umap_phase.pdf"),
       p_phase_umap, width = 8, height = 6)
ggsave(file.path(dir_plot, "04_cluster_phase_composition.pdf"),
       p_phase_comp, width = 10, height = 5)
message("[plot] 04_umap_phase.pdf + 04_cluster_phase_composition.pdf")

# ---- (B) Unbiased DE -------------------------------------------------------
future::plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)
message("[DE] FindAllMarkers (wilcox, positive only) ...")
markers <- FindAllMarkers(
  obj, assay = opt$assay, layer = "data", only.pos = TRUE,
  min.pct = opt$min_pct, logfc.threshold = opt$logfc, verbose = FALSE
)
if (nrow(markers) == 0)
  stop("[fatal] FindAllMarkers returned empty -- inspect [check] above.")
write_csv(markers, file.path(opt$outdir, "04_markers_all.csv"))

top_tbl <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = opt$top_n, with_ties = FALSE) |>
  ungroup()
write_csv(top_tbl, file.path(opt$outdir, "04_markers_top.csv"))
message(glue("[DE] {nrow(markers)} rows; top {opt$top_n}/cluster saved"))

dd_feats <- top_tbl |> distinct(gene) |> pull(gene)
p_dd <- DotPlot(obj, features = dd_feats) + RotatedAxis() +
  labs(title = glue("Top {opt$top_n} data-driven markers / cluster ({res_col})")) +
  theme(axis.text.x = element_text(size = 7))
ggsave(file.path(dir_plot, "04_dotplot_datadriven.pdf"),
       p_dd, width = max(12, 0.18 * length(dd_feats)), height = 6)
message("[plot] 04_dotplot_datadriven.pdf")

# ---- (C) Canonical panel dotplot -------------------------------------------
p_canon <- DotPlot(obj, features = panel, cluster.idents = FALSE) +
  RotatedAxis() +
  labs(title = glue("Canonical lineage panel ({res_col})")) +
  theme(axis.text.x = element_text(size = 7),
        strip.text.x  = element_text(size = 7, angle = 90))
ggsave(file.path(dir_plot, "04_dotplot_canonical.pdf"),
       p_canon, width = 22, height = 6.5)
message("[plot] 04_dotplot_canonical.pdf")

# ---- (C) Lineage module scores ---------------------------------------------
score_cols <- character(0)
for (lin in names(panel)) {
  obj <- AddModuleScore(obj, features = list(panel[[lin]]),
                        name = glue("{lin}_"), ctrl = 50, seed = 42)
  score_cols[lin] <- glue("{lin}_1")
}
score_mat <- obj@meta.data |>
  as_tibble() |>
  select(cluster = all_of(res_col), all_of(unname(score_cols))) |>
  group_by(cluster) |>
  summarise(across(everything(), mean), .groups = "drop")
colnames(score_mat) <- c("cluster", names(score_cols))
write_csv(score_mat, file.path(opt$outdir, "04_lineage_scores_by_cluster.csv"))

long <- score_mat |>
  pivot_longer(-cluster, names_to = "lineage", values_to = "score") |>
  group_by(lineage) |>
  mutate(z = as.numeric(scale(score))) |>
  ungroup() |>
  mutate(cluster = factor(cluster, levels = levels(Idents(obj))),
         lineage  = factor(lineage,  levels = names(panel)))
best <- long |> group_by(cluster) |> slice_max(z, n = 1) |> ungroup()
message("[first-pass] best-scoring lineage per cluster:")
print(best |> select(cluster, lineage, z) |> arrange(cluster), n = n_clust)

p_heat <- ggplot(long, aes(lineage, cluster, fill = z)) +
  geom_tile(color = "grey90") +
  geom_text(data = best, aes(label = "*"), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0) +
  labs(x = NULL, y = glue("cluster ({res_col})"), fill = "z",
       title = "Lineage module scores by cluster  (cell-cycle removed)",
       subtitle = "* = top lineage per cluster") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(dir_plot, "04_lineage_score_heatmap.pdf"),
       p_heat, width = 12, height = 7)
message("[plot] 04_lineage_score_heatmap.pdf")

saveRDS(obj, file.path(opt$outdir, "objects", "04_scored.rds"))
message("[done] Step 2 complete. Check [panel] dropped lines above, then review ",
        "04_umap_phase.pdf before interpreting cluster identities -- high S/G2M ",
        "clusters may be cycling progenitor states rather than distinct lineages.")