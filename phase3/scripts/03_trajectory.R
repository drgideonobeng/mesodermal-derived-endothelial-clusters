#!/usr/bin/env Rscript
# 03_trajectory.R
# Clean developmental trajectory of the endothelial compartment: progression from
# the Angioblast progenitor to the mature endothelial fates, with the branching
# structure and the connectivity between beds made explicit. No SHF / tip / origin
# logic here -- this script answers only "what is the lineage topology?"
#
# ── THE FIX (read before running) ──────────────────────────────────────────────
# The previous version of this script gave a biologically implausible topology
# (arterial-from-venous: Angioblast -> RA_producing -> Venous -> Posterior_
# sinusoidal, THEN splitting to Arterial_PAA). Root cause: it fed Slingshot bed
# labels BORROWED from the integrated (fastMNN) object, joined onto the
# UNCORRECTED geometry by barcode. Slingshot needs cluster labels that are
# coherent in the embedding it actually fits on -- labels from one embedding
# bolted onto another's geometry is a mismatch, not a trajectory.
#
# THE FIX: labels now come NATIVELY from this object's own clustering
# (`subcluster`, produced by 01 with --integrate FALSE), mapped to identity via
# the canonical-marker annotation already verified in 03a_traj_native_markers.R.
# No cross-embedding join, no --annot input -- the integrated object is not used
# here at all. This freezes a SECOND map (native_bed_map), independent of 02's
# integrated bed_map; the two use different cluster numberings and are not
# interchangeable.
#
# ── REMAINING CAVEATS (the fix above does not resolve these) ───────────────────
# • Stage gaps: the uncorrected embedding still has stage gaps (E8.0 | E8.25 |
#   late-merged), so a bed can fragment across UMAP islands and Slingshot
#   bridges the gaps. Branch TOPOLOGY (which beds connect, where it splits) is
#   informative; pseudotime VALUES across the early gap are interpolated, not
#   sampled. Read within-late-island ordering as reliable, cross-gap ordering
#   as approximate.
# • Ambient Hb in cluster 0: native cluster 0 (Venous_capillary_E105) carries
#   severe ambient hemoglobin contamination (Hba/Hbb 68-86% positive). 01's HVG
#   selection did NOT exclude Hb genes, so cluster 0's discreteness as a native
#   cluster could be partly an Hb-driven artifact rather than real biology. If
#   the topology below still looks implausible around cluster 0 specifically,
#   that -- not label-geometry mismatch -- is the next suspect. --drop_native
#   lets you refit excluding it as a quick diagnostic; the durable fix would be
#   re-running 01 (--integrate FALSE) with Hb genes excluded from HVGs first.
#
# Input : results/phase3/<lineage>/uncorrected/01_endo_subset.rds  (native `subcluster`)
# Output: 03_trajectory_lineages.pdf | 03_trajectory_pseudotime.pdf |
#         03_native_bed_umap.pdf     | 03_lineage_structure.txt |
#         03_trajectory_object.rds

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(tidyverse); library(glue); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage",     type = "character", default = "endothelium",
              help = "phase-3 namespace under results/phase3/<lineage>/"),
  make_option("--traj",        type = "character", default = NULL,
              help = "uncorrected subset; default results/phase3/<lineage>/uncorrected/01_endo_subset.rds"),
  make_option("--outdir",      type = "character", default = NULL,
              help = "default results/phase3/<lineage>/uncorrected"),
  make_option("--root",        type = "character", default = "Angioblast_E80"),
  make_option("--n_dims",      type = "integer",   default = 20L),
  make_option("--drop_native", type = "character", default = "",
              help = "comma-separated native subcluster ID(s) to exclude before fitting (diagnostic; e.g. '0' for the Hb-contaminated cluster)")
)))
set.seed(42)

# ── Resolve paths from --lineage ────────────────────────────────────────────────
base <- file.path("results/phase3", opt$lineage)
if (is.null(opt$traj))   opt$traj   <- file.path(base, "uncorrected", "01_endo_subset.rds")
if (is.null(opt$outdir)) opt$outdir <- file.path(base, "uncorrected")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

traj <- readRDS(opt$traj)
stopifnot("subcluster" %in% colnames(traj@meta.data))

# ── Freeze the NATIVE bed annotation (independent of 02's integrated bed_map) ─
# Identity by canonical markers only (integrin excluded), verified in
# 03a_traj_native_markers.R. Valid ONLY for this clustering (01, --integrate
# FALSE, locked seed/resolution) -- the guard below fails loudly if it differs.
native_bed_map <- c(
  "0" = "Venous_capillary_E105",   # Hba/Hbb 68-86% positive -- ambient contamination, see caveat above
  "1" = "Posterior_sinusoidal_E80",
  "2" = "Venous_E80",
  "3" = "Progenitor_E825",
  "4" = "RA_producing_E80",
  "5" = "Arterial_E95",            # Mecom, Bmx
  "6" = "Angioblast_E80",          # root
  "7" = "Arterial_E825"            # Tmem100
)
miss <- setdiff(levels(traj$subcluster), names(native_bed_map))
if (length(miss) > 0)
  stop(glue("[freeze] native subcluster(s) not mapped: {paste(miss, collapse=', ')} -- ",
            "the clustering in --traj differs from the one this map was built for. Re-verify via 03a."))
traj$bed_native <- factor(unname(native_bed_map[as.character(traj$subcluster)]),
                          levels = unname(native_bed_map))

if (!opt$root %in% levels(traj$bed_native))
  stop(glue("[traj] root '{opt$root}' not a native bed label; have: ",
            "{paste(levels(traj$bed_native), collapse = ', ')}"))

# ── Optional diagnostic drop (e.g. the Hb-contaminated cluster) ───────────────
drop_ids <- str_split(opt$drop_native, ",")[[1]] |> str_trim() |> discard(~ .x == "")
if (length(drop_ids) > 0) {
  drop_names <- native_bed_map[drop_ids]
  message(glue("[drop] excluding native cluster(s) {paste(drop_ids, collapse=', ')} ",
               "({paste(drop_names, collapse=', ')}) before fitting"))
  traj <- subset(traj, subset = !subcluster %in% drop_ids)
  traj$bed_native <- droplevels(traj$bed_native)
}

# ── QC: native cluster identity is self-consistent with its own geometry ──────
# (sanity check that the fix is doing what it claims -- labels and embedding
# now come from the same object/clustering, unlike the prior version)
p_clu <- DimPlot(traj, reduction = "umap_sub", group.by = "subcluster",
                 label = TRUE, repel = TRUE) + ggtitle("Native subcluster (integer)") + NoLegend()
p_bed <- DimPlot(traj, reduction = "umap_sub", group.by = "bed_native",
                 label = TRUE, repel = TRUE, label.size = 2.8) +
  ggtitle("Native bed identity (canonical markers)") + NoLegend()
ggsave(out("03_native_bed_umap.pdf"), p_clu | p_bed, width = 14, height = 6)

# -- Slingshot: lineage inference + pseudotime on the uncorrected PCA -----------
if (!"pca_sub" %in% Reductions(traj))
  stop("[traj] pca_sub reduction missing -- was the uncorrected object built with --integrate FALSE?")
dims <- seq_len(min(opt$n_dims, ncol(Embeddings(traj, "pca_sub"))))
rd   <- Embeddings(traj, "pca_sub")[, dims]

sds <- slingshot(rd, clusterLabels = as.character(traj$bed_native), start.clus = opt$root)
pt  <- slingPseudotime(sds)                          # cells x lineages
colnames(pt) <- paste0("pt_L", seq_len(ncol(pt)))
traj$pseudotime_mean <- rowMeans(pt, na.rm = TRUE)   # collapsed maturation axis
for (i in seq_len(ncol(pt))) traj[[colnames(pt)[i]]] <- pt[, i]

lin <- slingLineages(sds)

# -- Report the branching structure (which beds connect, where it splits) ------
struct <- imap_chr(lin, ~ glue("L{.y}: {paste(.x, collapse = ' -> ')}"))
writeLines(c(glue("Root: {opt$root}"),
             glue("Lineages inferred: {length(lin)}"),
             glue("Terminal fates: {paste(map_chr(lin, ~ tail(.x, 1)), collapse = ', ')}"),
             "", struct),
           out("03_lineage_structure.txt"))
message("[slingshot] lineage structure:"); walk(struct, ~ message("  ", .x))
message("[check] expect the arterial leg as Angioblast_E80 -> Arterial_E825 -> Arterial_E95, ",
        "NOT routed through a venous bed -- if it still routes through venous, see the Hb caveat above.")

# -- Bed centroids + inferred branch edges in UMAP space (the connectivity tree) -
um <- Embeddings(traj, "umap_sub") |> as_tibble(.name_repair = ~ c("u1","u2")) |>
  mutate(bed = traj$bed_native, pseudotime = traj$pseudotime_mean)
cent  <- um |> group_by(bed) |>
  summarise(u1 = median(u1), u2 = median(u2), .groups = "drop")
edges <- map_dfr(lin, ~ tibble(from = head(.x, -1), to = tail(.x, -1))) |> distinct() |>
  left_join(cent, by = c("from" = "bed")) |> rename(x = u1, y = u2) |>
  left_join(cent, by = c("to"   = "bed")) |> rename(xend = u1, yend = u2)

# -- Plot 1: lineage tree over the bed UMAP (branches + connectivity) -----------
p_lin <- ggplot(um, aes(u1, u2)) +
  geom_point(aes(colour = bed), size = .5, alpha = .7) +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = .55, colour = "grey25",
               arrow = arrow(length = unit(.12, "cm"), type = "closed")) +
  geom_point(data = cent, aes(u1, u2), size = 3.2, shape = 21,
             fill = "white", colour = "black") +
  geom_text(data = cent, aes(u1, u2, label = bed), size = 2.6, vjust = -1.2) +
  labs(title = glue("Endothelial lineage structure (root: {opt$root})"),
       subtitle = "cells coloured by native bed; nodes = bed centroids; arrows = inferred branch tree",
       colour = NULL) +
  theme_minimal(base_size = 11)
ggsave(out("03_trajectory_lineages.pdf"), p_lin, width = 9, height = 7)

# -- Plot 2: pseudotime maturation gradient over the UMAP ----------------------
p_pt <- ggplot(um, aes(u1, u2, colour = pseudotime)) +
  geom_point(size = .6, alpha = .85) +
  scale_colour_viridis_c(na.value = "grey85") +
  labs(title = glue("Pseudotime: {opt$root} -> mature fates"),
       subtitle = "mean across lineages; cross-gap ordering interpolated, not sampled",
       colour = "pseudotime") +
  theme_minimal(base_size = 11)
ggsave(out("03_trajectory_pseudotime.pdf"), p_pt, width = 8, height = 6)

saveRDS(list(sds = sds, traj = traj, lineages = lin), out("03_trajectory_object.rds"))
message(glue("[done] -> {out('03_trajectory_object.rds')}"))