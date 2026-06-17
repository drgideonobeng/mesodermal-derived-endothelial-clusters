#!/usr/bin/env Rscript
# 06_annotate.R
# Canonical lineage panel: grouped dotplot + per-lineage module scores.
# Panel markers are grouped by lineage; absent genes are dropped and reported.
# Input : 04_cc_scored.rds  +  05_markers_top.csv
# Output: 06_dotplot_canonical.pdf | 06_score_heatmap.pdf | 06_module_scores.csv

# ── PANEL NOTE ────────────────────────────────────────────────────────────────
# Markers were compiled from canonical developmental literature. For publication,
# verify against the original Nomaru et al. 2021 (GSE158941) marker tables and
# cross-check with CellMarker 2.0 or PanglaoDB. Absent genes are silently
# dropped by Seurat — the [panel] log lines below list exactly which ones.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type = "character"),
  make_option("--markers",    type = "character"),
  make_option("--resolution", type = "double", default = 0.4)
)))

res_col <- glue("clust_res{opt$resolution}")

# ── Lineage panel (identity markers only — no cell-cycle genes) ───────────────
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
  Lymphatic_EC     = c("Prox1","Lyve1","Pdpn","Flt4"),   # expected near-null at E8-10.5
  Endocardium      = c("Npr3","Nfatc1"),
  Erythroid_blood  = c("Hba-x","Hbb-bh1","Hbb-y","Hba-a1","Gata1","Klf1","Runx1"),
  LPM_limb         = c("Hand2","Foxf1","Prrx1","Tbx4"),
  Epicardium       = c("Wt1","Tbx18","Upk3b"),
  Smooth_mural     = c("Acta2","Tagln","Myh11","Pdgfrb"),
  Mesenchyme_EMT   = c("Twist1","Snai1","Postn","Col1a1"),
  Notochord_axial  = c("T","Tbxt","Noto","Shh"),         # T and Tbxt: same gene, annotation varies
  Neural_ectoderm  = c("Sox2","Pax6","Tubb3","Sox1"),
  Endoderm         = c("Foxa2","Sox17","Epcam","Cldn6","Afp")
)

# ── Load + filter panel ───────────────────────────────────────────────────────
obj     <- readRDS(opt$input)
present <- rownames(obj)
Idents(obj) <- res_col

panel_raw <- panel
panel     <- map(panel, ~ .x[.x %in% present])
dropped   <- map2(panel_raw, panel, ~ setdiff(.x, .y))
dropped   <- keep(dropped, ~ length(.x) > 0)
if (length(dropped) > 0) {
  message("[panel] genes not in object (dropped):")
  iwalk(dropped, ~ message(glue("  {.y}: {paste(.x, collapse=', ')}")))
}
panel <- keep(panel, ~ length(.x) > 0)
message(glue("[panel] {length(unlist(panel))} genes / {length(panel)} lineages"))

# ── Canonical dotplot ─────────────────────────────────────────────────────────
p_dot <- DotPlot(obj, features = panel) + RotatedAxis() +
  labs(title = glue("Canonical lineage panel ({res_col})")) +
  theme(axis.text.x = element_text(size = 7),
        strip.text.x = element_text(size = 7, angle = 90))
ggsave("06_dotplot_canonical.pdf", p_dot, width = 22, height = 6.5)

# ── Module scores (z-score per lineage across clusters) ───────────────────────
for (lin in names(panel))
  obj <- AddModuleScore(obj, features = list(panel[[lin]]),
                        name = glue("{lin}_"), ctrl = 50, seed = 42)

score_cols <- setNames(glue("{names(panel)}_1"), names(panel))
score_mat  <- obj@meta.data |>
  as_tibble() |>
  select(cluster = all_of(res_col), all_of(unname(score_cols))) |>
  group_by(cluster) |>
  summarise(across(everything(), mean), .groups = "drop")
colnames(score_mat) <- c("cluster", names(panel))
write_csv(score_mat, "06_module_scores.csv")

long <- score_mat |>
  pivot_longer(-cluster, names_to = "lineage", values_to = "score") |>
  group_by(lineage) |>
  mutate(z       = as.numeric(scale(score)),
         cluster = factor(cluster, levels = levels(Idents(obj))),
         lineage = factor(lineage, levels = names(panel))) |>
  ungroup()

best <- long |> group_by(cluster) |> slice_max(z, n = 1) |> ungroup()
message("[annotate] first-pass best lineage per cluster:")
print(select(arrange(best, cluster), cluster, lineage, z), n = nlevels(Idents(obj)))

p_heat <- ggplot(long, aes(lineage, cluster, fill = z)) +
  geom_tile(color = "grey90") +
  geom_text(data = best, aes(label = "*"), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0) +
  labs(x = NULL, y = glue("cluster ({res_col})"), fill = "z",
       title  = "Lineage module scores by cluster",
       subtitle = "z = SD from mean across clusters | * = top lineage") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("06_score_heatmap.pdf", p_heat, width = 12, height = 7)

message("[annotate] done")
