#!/usr/bin/env Rscript
# 11_trajectory.R
# Clean developmental trajectory of the endothelial compartment: progression from
# the Angioblast progenitor to the mature endothelial fates, with the branching
# structure and the connectivity between beds made explicit. No SHF / tip / origin
# logic here -- this script answers only "what is the lineage topology?"
#
# Slingshot is fit on the UNCORRECTED PCA (pca_sub): fastMNN flattened the
# E8.0->E10.5 axis, so the uncorrected embedding is the correct substrate for
# developmental ordering. Lineage inference + pseudotime come from the PCA; the
# branch tree and maturation gradient are visualised on the UMAP.
#
# Bed labels are carried by barcode from the annotated (integrated) object, so
# identity comes from the robust clustering while geometry comes from the
# uncorrected embedding (same 1049 cells in both).
#
# CAVEAT: the uncorrected embedding has stage gaps (E8.0 | E8.25 | late-merged),
# so a bed can fragment across islands and Slingshot bridges the gaps. The branch
# TOPOLOGY (which beds connect, where it splits) is informative; pseudotime VALUES
# across the early gap are interpolated, not sampled. Read ordering WITHIN the
# continuous late island as reliable and cross-gap ordering as approximate. The
# late-window refit (next script) is what makes gene-vs-pseudotime quantitative.
#
# Input : --traj  results/phase2a/endo_traj/09_endo_subset.rds    (uncorrected)
#         --annot results/phase2a/endo_only/10_endo_annotated.rds  (bed labels)
# Output: 11_trajectory_lineages.pdf | 11_trajectory_pseudotime.pdf |
#         11_lineage_structure.txt   | 11_trajectory_object.rds

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(tidyverse); library(glue); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--traj",   type = "character",
              default = "results/phase2a/endo_traj/09_endo_subset.rds"),
  make_option("--annot",  type = "character",
              default = "results/phase2a/endo_only/10_endo_annotated.rds"),
  make_option("--outdir", type = "character", default = "results/phase2a/endo_traj"),
  make_option("--root",   type = "character", default = "Angioblast"),
  make_option("--n_dims", type = "integer",   default = 20L)
)))
set.seed(42)
out <- function(f) file.path(opt$outdir, f)

traj  <- readRDS(opt$traj)
annot <- readRDS(opt$annot)

# -- Carry bed identity over by barcode (same 1049 endothelial cells) ----------
if (!setequal(colnames(traj), colnames(annot)))
  stop("[traj] cell barcodes differ between the uncorrected and annotated objects")
traj$bed <- factor(annot$bed[match(colnames(traj), colnames(annot))],
                   levels = levels(annot$bed))
if (!opt$root %in% levels(traj$bed))
  stop(glue("[traj] root '{opt$root}' not a bed label; have: ",
            "{paste(levels(traj$bed), collapse = ', ')}"))

# -- Slingshot: lineage inference + pseudotime on the uncorrected PCA -----------
if (!"pca_sub" %in% Reductions(traj))
  stop("[traj] pca_sub reduction missing -- was the uncorrected object built with --integrate FALSE?")
dims <- seq_len(min(opt$n_dims, ncol(Embeddings(traj, "pca_sub"))))
rd   <- Embeddings(traj, "pca_sub")[, dims]

sds <- slingshot(rd, clusterLabels = as.character(traj$bed), start.clus = opt$root)
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
           out("11_lineage_structure.txt"))
message("[slingshot] lineage structure:"); walk(struct, ~ message("  ", .x))

# -- Bed centroids + inferred branch edges in UMAP space (the connectivity tree) -
um <- Embeddings(traj, "umap_sub") |> as_tibble(.name_repair = ~ c("u1","u2")) |>
  mutate(bed = traj$bed, pseudotime = traj$pseudotime_mean)
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
       subtitle = "cells coloured by bed; nodes = bed centroids; arrows = inferred branch tree",
       colour = NULL) +
  theme_minimal(base_size = 11)
ggsave(out("11_trajectory_lineages.pdf"), p_lin, width = 9, height = 7)

# -- Plot 2: pseudotime maturation gradient over the UMAP ----------------------
p_pt <- ggplot(um, aes(u1, u2, colour = pseudotime)) +
  geom_point(size = .6, alpha = .85) +
  scale_colour_viridis_c(na.value = "grey85") +
  labs(title = glue("Pseudotime: {opt$root} -> mature fates"),
       subtitle = "mean across lineages; cross-gap ordering interpolated, not sampled",
       colour = "pseudotime") +
  theme_minimal(base_size = 11)
ggsave(out("11_trajectory_pseudotime.pdf"), p_pt, width = 8, height = 6)

saveRDS(list(sds = sds, traj = traj, lineages = lin), out("11_trajectory_object.rds"))
message(glue("[done] -> {out('11_trajectory_object.rds')}"))