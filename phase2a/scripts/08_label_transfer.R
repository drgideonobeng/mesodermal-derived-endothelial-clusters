#!/usr/bin/env Rscript
# 08_label_transfer.R
# Independent identity call via reference label transfer (Pijuan-Sala 2019 mouse
# gastrulation atlas) onto the integrated WT manifold. Primary purpose: resolve
# the two clusters the canonical panel could not (Anterior_unresolved=3,
# FHF_unresolved=14) and confirm the frozen lineages against an external atlas.
#
# REFERENCE: MouseGastrulationData::EmbryoAtlasData() spans E6.5-E8.5, matching
# the EARLY end of this timecourse — exactly where the unresolved clusters sit
# (both E8.0/E8.25). Late clusters (E9.5/E10.5) are OUT OF REFERENCE RANGE, so
# treat their transferred labels with caution; a MOCA pass would complement.
#
# DEPENDENCY: MouseGastrulationData (BiocManager::install("MouseGastrulationData"))
# and a one-time ExperimentHub download (cached afterwards; needs internet).
# Run standalone first to do the download + inspect the crosstab before wiring
# into the pipeline:
#   Rscript scripts/08_label_transfer.R --input results/phase2a/objects/07_annotated.rds
#
# Input : 07_annotated.rds
# Output: 08_label_transfer_object.rds | 08_transfer_crosstab.csv |
#         08_transfer_heatmap.pdf | 08_predicted_umap.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(MouseGastrulationData)
  library(tidyverse)
  library(glue)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input",  type = "character"),
  make_option("--n_dims", type = "integer",  default = 30L),
  make_option("--stages", type = "character", default = "E8.0,E8.25,E8.5",
              help = "reference stages to include (comma-separated)")
)))

stages <- str_split(opt$stages, ",")[[1]] |> str_trim()

# ── Reference: Pijuan-Sala, early stages, real cells only ─────────────────────
message(glue("[transfer] loading EmbryoAtlasData (stages: {paste(stages, collapse=', ')}) ..."))
ref_sce <- EmbryoAtlasData()
ref_sce <- ref_sce[, ref_sce$stage %in% stages &
                     !ref_sce$doublet & !ref_sce$stripped &
                     !is.na(ref_sce$celltype)]

# rownames are Ensembl IDs -> map to gene symbols to match the query
rownames(ref_sce) <- rowData(ref_sce)$SYMBOL
ref_sce <- ref_sce[!is.na(rownames(ref_sce)) & !duplicated(rownames(ref_sce)), ]

ref <- CreateSeuratObject(counts    = counts(ref_sce),
                          meta.data = as.data.frame(colData(ref_sce))) |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
message(glue("[transfer] reference: {ncol(ref)} cells / {n_distinct(ref$celltype)} cell types"))

# ── Query (frozen object), log-normalised on RNA ──────────────────────────────
query <- readRDS(opt$input)
DefaultAssay(query) <- "RNA"
query[["RNA"]] <- JoinLayers(query[["RNA"]])
query <- NormalizeData(query, verbose = FALSE)

# ── Anchor-based transfer (expression space, not the MNN embedding) ───────────
anchors <- FindTransferAnchors(reference = ref, query = query,
                               features = VariableFeatures(ref),
                               dims = 1:opt$n_dims, reduction = "pcaproject",
                               verbose = FALSE)
pred <- TransferData(anchorset = anchors, refdata = ref$celltype,
                     dims = 1:opt$n_dims, verbose = FALSE)
query$predicted_celltype <- pred$predicted.id
query$predicted_score    <- pred$prediction.score.max

# ── Crosstab: frozen lineage vs transferred reference label ───────────────────
cross <- query@meta.data |>
  as_tibble() |>
  count(lineage, predicted_celltype) |>
  group_by(lineage) |>
  mutate(frac = n / sum(n)) |>
  ungroup()
write_csv(cross, "08_transfer_crosstab.csv")

top_pred <- cross |> group_by(lineage) |> slice_max(frac, n = 1, with_ties = FALSE) |> ungroup()
message("[transfer] dominant reference label per frozen lineage:")
print(arrange(dplyr::select(top_pred, lineage, predicted_celltype, frac), lineage),
      n = nrow(top_pred))

# ── Plots ─────────────────────────────────────────────────────────────────────
p_heat <- ggplot(cross, aes(predicted_celltype, lineage, fill = frac)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  labs(title = "Frozen lineage vs Pijuan-Sala transferred label",
       x = "predicted reference cell type", y = "frozen lineage", fill = "frac") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave("08_transfer_heatmap.pdf", p_heat, width = 14, height = 7)

p_umap <- DimPlot(query, reduction = "umap.mnn", group.by = "predicted_celltype",
                  label = TRUE, repel = TRUE, label.size = 2.5) +
  ggtitle("Transferred reference labels (Pijuan-Sala)") + NoLegend()
ggsave("08_predicted_umap.pdf", p_umap, width = 11, height = 8)

saveRDS(query, "08_label_transfer_object.rds")
message("[transfer] done → 08_label_transfer_object.rds")