#!/usr/bin/env Rscript
# 02_shf_modules.R
# Harden the non-SHF / compensatory-cell leads from the Arterial_PAA (subcluster 3)
# SHF split. The single-gene DE (Esm1/Adm/Hoxb up in SHF-low) did NOT survive
# multiple-testing correction — underpowered with only ~56 SHF-low cells. Here we
# aggregate the signal into MODULE SCORES and test the *programs*, which has far
# more power than any individual gene.
#
# ── ANTI-DOUBLE-DIPPING (read this) ───────────────────────────────────────────
# The hypothesis ("non-SHF cells are tip-like / more posterior") came FROM the DE
# on this same split. Testing a score that re-uses those exact genes on the same
# cells is circular and inflates significance. So both programs are built from
# CANONICAL / literature genes and DELIBERATELY EXCLUDE the DE hits:
#   Tip program     — excludes Esm1, Adm   (canonical sprouting markers instead)
#   PostHox program — excludes Hoxb2, Hoxb5 (broad arch-level paralog set instead)
# A positive result is then genuine corroboration, not the same finding restated.
#
# Also freezes the endothelial BED annotation (subcluster int -> bed name). The
# map is valid ONLY for the frozen clustering in --input; a guard fails loudly if
# the clustering differs.
#
# Input : results/phase3/endothelium/integrated/01_endo_subset.rds  (integrated, 6 beds)
# Output: 02_endo_annotated.rds | 02_module_violins.pdf | 02_module_tests.csv

suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(glue); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lineage",   type = "character", default = "endothelium",
              help = "phase-3 namespace under results/phase3/<lineage>/"),
  make_option("--input",     type = "character", default = NULL,
              help = "integrated subset object; default results/phase3/<lineage>/integrated/01_endo_subset.rds"),
  make_option("--outdir",    type = "character", default = NULL,
              help = "default results/phase3/<lineage>/integrated"),
  make_option("--focus_bed", type = "character", default = "Arterial_PAA",
              help = "bed within which to test SHF-high vs SHF-low")
)))
set.seed(42)

# ── Resolve paths from --lineage (override with --input / --outdir) ────────────
base <- file.path("results/phase3", opt$lineage)
if (is.null(opt$input))  opt$input  <- file.path(base, "integrated", "01_endo_subset.rds")
if (is.null(opt$outdir)) opt$outdir <- file.path(base, "integrated")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

sub <- readRDS(opt$input)
stopifnot("subcluster" %in% colnames(sub@meta.data))

# ── Freeze the endothelial bed annotation (valid for THIS clustering only) ────
bed_map <- c("0" = "Venous", "1" = "Posterior_sinusoidal", "2" = "RA_producing",
             "3" = "Arterial_PAA", "4" = "Angioblast", "5" = "Venous_lymphatic")
miss <- setdiff(levels(sub$subcluster), names(bed_map))
if (length(miss))
  stop(glue("[freeze] subcluster(s) not mapped: {paste(miss, collapse=', ')} — ",
            "the clustering in --input differs from the annotated one. Re-check the object."))
sub$bed <- factor(unname(bed_map[as.character(sub$subcluster)]),
                  levels = unname(bed_map))

# ── Module-score gene sets ────────────────────────────────────────────────────
shf_genes <- c("Isl1","Tbx1","Mef2c","Fgf8","Fgf10","Six2","Hand2")        # lineage definers
tip_genes <- c("Apln","Angpt2","Pgf","Kcne3","Igfbp3","Cxcr4","Dll4",       # canonical tip/sprouting
               "Cd34","Nid2","Lxn")                                          #   (Esm1/Adm EXCLUDED)
posthox   <- c("Hoxa1","Hoxa2","Hoxa3","Hoxa4","Hoxb1","Hoxb3","Hoxb4",     # arch-level Hox
               "Hoxb6","Hoxb8")                                              #   (Hoxb2/Hoxb5 EXCLUDED)

gene_sets <- list(SHF = shf_genes, Tip = tip_genes, PostHox = posthox)
gene_sets <- lapply(gene_sets, \(g) g[g %in% rownames(sub)])
walk2(names(gene_sets), gene_sets,
      ~ message(glue("[score] {.x}: {length(.y)} genes present ({paste(.y, collapse=', ')})")))

sub <- AddModuleScore(sub, features = gene_sets, name = "mod", ctrl = 50, seed = 42)
sub$SHF_score     <- sub$mod1   # list order: SHF, Tip, PostHox
sub$Tip_score     <- sub$mod2
sub$PostHox_score <- sub$mod3
sub$mod1 <- sub$mod2 <- sub$mod3 <- NULL

# ── SHF status (>0 = above expression-matched background) ─────────────────────
sub$shf_status <- ifelse(sub$SHF_score > 0, "SHF_high", "SHF_low")

foc_cells <- colnames(sub)[sub$bed == opt$focus_bed]
md <- sub@meta.data[foc_cells, ] |> as_tibble()
message(glue("[focus] {opt$focus_bed}: {sum(md$shf_status=='SHF_high')} SHF_high / ",
             "{sum(md$shf_status=='SHF_low')} SHF_low"))

# ── Test each program across the SHF split (Wilcoxon + Cliff's delta) ─────────
# Cliff's delta in [-1,1]; positive => SHF_low scores higher (our directional
# expectation for both tip and posterior-Hox).
cliffs <- function(lo, hi) mean(outer(lo, hi, ">")) - mean(outer(lo, hi, "<"))
test_one <- function(score) {
  hi <- md[[score]][md$shf_status == "SHF_high"]
  lo <- md[[score]][md$shf_status == "SHF_low"]
  w  <- wilcox.test(lo, hi)
  tibble(program = score, median_hi = median(hi), median_lo = median(lo),
         delta_lo_minus_hi = median(lo) - median(hi),
         cliffs_delta = cliffs(lo, hi), p = w$p.value)
}
res <- map_df(c("Tip_score","PostHox_score"), test_one) |>
  mutate(p_adj = p.adjust(p, "BH"))
print(res); write_csv(res, out("02_module_tests.csv"))

# ── Visualise the split ───────────────────────────────────────────────────────
vln <- md |> select(shf_status, Tip_score, PostHox_score) |>
  pivot_longer(-shf_status, names_to = "program", values_to = "score") |>
  ggplot(aes(shf_status, score, fill = shf_status)) +
  geom_violin(scale = "width", alpha = .7) +
  geom_boxplot(width = .15, outlier.size = .4) +
  facet_wrap(~ program, scales = "free_y") +
  labs(title = glue("Independent programs across SHF split in {opt$focus_bed}"),
       subtitle = "Tip/PostHox genes exclude the DE hits — independent test of the lead") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
ggsave(out("02_module_violins.pdf"), vln, width = 8, height = 4)

saveRDS(sub, out("02_endo_annotated.rds"))
message(glue("[done] bed + scores frozen -> {out('02_endo_annotated.rds')}"))
