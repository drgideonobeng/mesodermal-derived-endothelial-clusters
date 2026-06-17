#!/usr/bin/env Rscript
# 04_cellcycle.R
# Score cell cycle phase (Tirosh S/G2M gene sets, mouse-converted) and
# produce phase UMAP + per-cluster composition bar. Scoring is always run;
# regression is deliberately omitted (see note below).
# Input : 03_clustered.rds
# Output: 04_cc_scored.rds | 04_umap_phase.pdf |
#         04_phase_composition.pdf | 04_phase_by_cluster.csv

# ── CELL CYCLE REGRESSION NOTE ───────────────────────────────────────────────
# Regression (vars.to.regress = c("S.Score","G2M.Score") in ScaleData) removes
# cell-cycle variation before embedding and clustering. Use it when:
#   - Cell cycle is the dominant transcriptional axis masking biology of interest
#     (common in tumour or in-vitro datasets with mixed cycling/quiescent cells).
# Do NOT regress in developmental data like this timecourse:
#   - Rapidly dividing progenitors SHOULD differ from post-mitotic derivatives.
#     Regression collapses progenitor and differentiated states, erasing signal.
# fastMNN note: ScaleData regression affects PCA but NOT the MNN embedding
#   (which takes log-normalised counts directly). To remove cell-cycle influence
#   from the MNN space, exclude cc genes from the HVG list instead.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",      type = "character"),
  make_option("--resolution", type = "double", default = 0.4)
)))

res_col    <- glue("clust_res{opt$resolution}")
phase_cols <- c(G1 = "#4575b4", S = "#d73027", G2M = "#fdae61")

obj <- readRDS(opt$input)

# ── Score (human gene sets → mouse via str_to_title) ─────────────────────────
s_genes   <- str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
message(glue("[cellcycle] S: {sum(s_genes %in% rownames(obj))}/{length(s_genes)} genes found"))
message(glue("[cellcycle] G2M: {sum(g2m_genes %in% rownames(obj))}/{length(g2m_genes)} genes found"))

obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                        set.ident = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_umap <- DimPlot(obj, reduction = "umap.mnn", group.by = "Phase",
                  cols = phase_cols) +
  ggtitle("Cell cycle phase")

phase_comp <- obj@meta.data |>
  as_tibble() |>
  count(cluster = .data[[res_col]], Phase) |>
  group_by(cluster) |>
  mutate(frac = n / sum(n))

p_bar <- ggplot(phase_comp, aes(cluster, frac, fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = phase_cols) +
  labs(x = glue("cluster ({res_col})"), y = "fraction",
       title = "Cell cycle composition by cluster") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

ggsave("04_umap_phase.pdf",        p_umap, width = 8,  height = 6)
ggsave("04_phase_composition.pdf", p_bar,  width = 10, height = 5)
write_csv(phase_comp, "04_phase_by_cluster.csv")

saveRDS(obj, "04_cc_scored.rds")
message("[cellcycle] done → 04_cc_scored.rds")
