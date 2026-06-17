#!/usr/bin/env Rscript
# 07_annotate_freeze.R
# Freeze biological identities into a stable `lineage` metadata column so every
# downstream step references names, not Seurat's size-ranked cluster integers
# (which reshuffle whenever cell counts shift). The map below is the SINGLE
# pinned cluster->identity binding, from module scores + FindAllMarkers + the
# canonical dotplot at res 0.4 + reference label transfer. Clusters 3 and 14 
# particulary RESOLVED via reference label transfer (Script:08).
#
# NOTE: this map is valid for the CURRENT clustering. Re-running with the same
# inputs/seed is deterministic (numbering is stable), but if QC, resolution, or
# the upstream singlets change, re-verify the map — the guard below fails loudly
# if a cluster appears that isn't mapped.
#
# Input : 04_cc_scored.rds
# Output: 07_annotated.rds | 07_lineage_map.csv | 07_umap_lineage.pdf

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

res_col <- glue("clust_res{opt$resolution}")

# ── Cluster -> lineage map (res 0.4) ──────────────────────────────────────────
# Identities from convergent evidence: canonical markers + module scores,
# unbiased DE (FindAllMarkers), reference label transfer (08; Pijuan-Sala
# E8.0-8.5), and single-cell co-expression + developmental-stage logic for
# contested calls (NMP, primitive vs definitive erythroid, endothelium vs
# endoderm). Caveat: reference is E8.0-8.5, so late clusters (endocardium,
# epicardium, definitive blood) rest on markers, not transfer.
lineage_map <- c(
  "0"  = "Lateral_plate_mesoderm",          # Osr1, Foxf1, Fendrr
  "1"  = "Second_heart_field",              # Sfrp5, Mab21l2 (softer call)
  "2"  = "Cranial_mesenchyme",              # Alx3 (softer call)
  "3"  = "Cranial_paraxial_mesoderm",       # Otx2, Zic2 (anterior); ref Paraxial (08)
  "4"  = "Branchiomeric_muscle_progenitor", # Myf5+ in cardiopharyngeal field(Tbx1/Isl1/Pitx2+head bHLHs Msc/Tcf21; Pax3-independent; Myog low = progenitor
  "5"  = "Somitic_mesoderm",                # Meox1 93%, Pax3 48%; Tbx6/Msgn1 ~0 → not PSM
  "6"  = "Allantois",                       # Tbx4 84%, posterior Hox; ref Allantois (08); NMP ruled out
  "7"  = "Caudal_mesoderm",                 # Cdx1, Cdx4, Lhx1; NMP ruled out (Sox2/Bra co-exp ~noise)
  "8"  = "Endothelium",                     # Cldn5, Pecam1, Etv2 (Cldn6 cross-reactivity excluded)
  "9"  = "Primitive_erythroid",             # Hemgn, Hbb-y (embryonic globin) + E≤10.5 + Mesp1 yolk-sac origin
  "10" = "Pharyngeal_mesenchyme",           # Pax9, Sostdc1
  "11" = "Hematopoietic",                   # Rac2, Fermt3, F10
  "12" = "Cardiomyocyte",                   # Nppa, Myh6, Trdn
  "13" = "Endocardium",                     # Klf2, Bmx
  "14" = "Mesenchyme",                      # atypical DE; ref Mesenchyme 0.995 (08) — resolved
  "15" = "Epicardium"                       # Upk3b, Aldh1a1
)

obj      <- readRDS(opt$input)
clusters <- as.character(obj[[res_col]][, 1])

missing <- setdiff(unique(clusters), names(lineage_map))
if (length(missing) > 0)
  stop(glue("[freeze] clusters with no mapping: {paste(missing, collapse=', ')} — update lineage_map"))

obj$lineage <- factor(unname(lineage_map[clusters]), levels = unique(unname(lineage_map)))
Idents(obj) <- "lineage"
message(glue("[freeze] {nlevels(obj$lineage)} lineages across {ncol(obj)} cells"))

# ── Record the map + a labelled UMAP ──────────────────────────────────────────
tibble(cluster = names(lineage_map), lineage = unname(lineage_map)) |>
  write_csv("07_lineage_map.csv")

p <- DimPlot(obj, reduction = "umap.mnn", group.by = "lineage",
             label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle(glue("Frozen lineage identities (res {opt$resolution})")) + NoLegend()
ggsave("07_umap_lineage.pdf", p, width = 10, height = 8)

saveRDS(obj, "07_annotated.rds")
message("[freeze] done → 07_annotated.rds")
