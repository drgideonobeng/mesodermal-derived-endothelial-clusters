library(Seurat)
install.packages("SingleCellExperiment")
library(SingleCellExperiment)
library(tidyverse)
obj <- readRDS("~/bioinformatics/scRNAseq/paa-dev-trajectory/results/phase2a/objects/07_annotated.rds")
bra <- intersect(c("Tbxt","T"), rownames(obj))[1]      # which symbol your data uses
message("Brachyury symbol: ", bra)
d <- FetchData(obj, vars = c("Sox2", bra, "clust_res0.4"))
colnames(d) <- c("Sox2","Bra","cluster")
d |>
  filter(cluster %in% c("6","7")) |>
  mutate(state = case_when(
    Sox2 > 0 & Bra > 0 ~ "Sox2+ Bra+  (NMP)",
    Sox2 > 0           ~ "Sox2+ only",
    Bra  > 0           ~ "Bra+ only  (mesoderm)",
    TRUE               ~ "double-neg")) |>
  count(cluster, state) |>
  group_by(cluster) |> mutate(frac = n / sum(n)) |> arrange(cluster, desc(frac))

FeaturePlot(obj, c("Sox2", bra), blend = TRUE, reduction = "umap.mnn")

FetchData(obj, c("Tbx4", "clust_res0.4")) |>
  group_by(cluster = `clust_res0.4`) |>
  summarise(pct_Tbx4 = mean(Tbx4 > 0), mean_Tbx4 = mean(Tbx4)) |>
  arrange(desc(pct_Tbx4))

FetchData(obj, c("Hbb-bs", "clust_res0.4")) |>
  group_by(cluster = `clust_res0.4`) |>
  summarise(pct_Hbb_bs  = mean(`Hbb-bs` > 0),
            mean_Hbb_bs = mean(`Hbb-bs`)) |>
  arrange(desc(pct_Hbb_bs))

FetchData(obj, c("Meox1","Pax3","Tbx6","Msgn1","clust_res0.4")) |>
  group_by(cluster = `clust_res0.4`) |>
  summarise(across(c(Meox1,Pax3,Tbx6,Msgn1), ~ mean(.x > 0))) |>
  arrange(desc(Meox1))

genes <- c("Myf5","Myod1","Myog",                 # myogenic determination (which it is)
           "Tbx1","Tcf21","Msc","Isl1","Pitx2",   # branchiomeric / cardiopharyngeal
           "Pax3","Pax7","Lbx1")                   # somitic / trunk
FetchData(obj, c(genes, "clust_res0.4")) |>
  group_by(cluster = `clust_res0.4`) |>
  summarise(across(all_of(genes), ~ mean(.x > 0))) |>
  arrange(desc(Myf5))

obj <- readRDS("07_annotated.rds")
table(obj$lineage)   # 16 named levels, no NA, counts match the cluster sizes

library(Seurat); library(tidyverse)
t8 <- readRDS("~/bioinformatics/scRNAseq/paa-dev-trajectory/results/phase2a/objects/08_label_transfer_object.rds")  # predictions intact
o7 <- readRDS("~/bioinformatics/scRNAseq/paa-dev-trajectory/results/phase2a/objects/07_annotated.rds")              # corrected lineage

t8$lineage <- o7@meta.data[colnames(t8), "lineage"]                    # barcode-matched swap

cross <- t8@meta.data |> as_tibble() |>
  count(lineage, predicted_celltype) |>
  group_by(lineage) |> mutate(frac = n / sum(n)) |> ungroup()
write_csv(cross, "08_transfer_crosstab.csv")

library(glue)
library(fs)
root_path <- path("~", "bioinformatics", "scRNAseq", "paa-dev-trajectory")
setwd("~/bioinformatics/scRNAseq/paa-dev-trajectory")
file_path <- path("results", "phase2a", "endo_only", "09_endo_subset.rds")
sub <- readRDS(file_path)
nlevels(sub$subcluster); table(sub$subcluster)

DotPlot(sub, group.by = "subcluster",
        features = c("Isl1","Tbx1","Mef2c","Fgf8","Fgf10","Six2","Hand2",  # SHF / AHF program
                     "Cxcr4","Flt4","Plxnd1","Gbx2")) +                      # pharyngeal / arch-EC
  RotatedAxis()

# SHF / AHF lineage program — your Cre-line definers + core AHF
shf <- list(SHF = c("Isl1","Tbx1","Mef2c","Fgf8","Fgf10","Six2","Hand2"))
shf$SHF <- shf$SHF[shf$SHF %in% rownames(sub)]                 # keep present genes

sub <- AddModuleScore(sub, features = shf, name = "SHF", ctrl = 50, seed = 42)  # -> SHF1

# per-subcluster SHF score (which beds are SHF-derived?)
sub@meta.data |> as_tibble() |>
  group_by(subcluster) |>
  summarise(mean_SHF = mean(SHF1), pct_pos = mean(SHF1 > 0), n = n()) |>
  arrange(desc(mean_SHF)) |> print()

library(patchwork)
(VlnPlot(sub, "SHF1", group.by = "subcluster", pt.size = 0) + NoLegend()) |
  FeaturePlot(sub, "SHF1", reduction = "umap_sub")

art <- subset(sub, subcluster == "3")                          # arterial/PAA bed
art$shf_status <- ifelse(art$SHF1 > 0, "SHF_high", "SHF_low")  # >0 = above background
print(table(art$shf_status))

Idents(art) <- "shf_status"
de <- FindMarkers(art, ident.1 = "SHF_high", ident.2 = "SHF_low",
                  logfc.threshold = 0.25, min.pct = 0.1)
head(de, 30)

sub <- readRDS("results/phase2a/endo_only/10_endo_annotated.rds")
tip2 <- list(Tip = c("Apln","Angpt2","Pgf","Kcne3","Igfbp3","Cd34","Nid2","Lxn"))
tip2$Tip <- tip2$Tip[tip2$Tip %in% rownames(sub)]
sub <- AddModuleScore(sub, features = tip2, name = "tipc", ctrl = 50, seed = 42)
f <- sub@meta.data[sub$bed == "Arterial_PAA", ]
wilcox.test(f$tipc1[f$shf_status=="SHF_low"], f$tipc1[f$shf_status=="SHF_high"])

sub <- readRDS("results/phase2a/endo_only/10_endo_annotated.rds")
library(tidyverse)

ven <- list(Venous = c("Nr2f2","Nr2f1","Nrp2","Aplnr","Ephb4","Nt5e","Ackr1"))
ven$Venous <- ven$Venous[ven$Venous %in% rownames(sub)]
sub <- AddModuleScore(sub, features = ven, name = "Ven", ctrl = 50, seed = 42)
sub$Venous_score <- sub$Ven1; sub$Ven1 <- NULL

cliffs <- function(lo, hi) mean(outer(lo, hi, ">")) - mean(outer(lo, hi, "<"))

md <- FetchData(sub, c("Nr2f2","bed","shf_status","Venous_score")) |>
  filter(bed == "Arterial_PAA")
lo <- md$Venous_score[md$shf_status=="SHF_low"]; hi <- md$Venous_score[md$shf_status=="SHF_high"]
cat(sprintf("Venous score | low %.3f vs high %.3f | cliffs %.3f | p %.3g\n",
            median(lo), median(hi), cliffs(lo,hi), wilcox.test(lo,hi)$p.value))
nlo <- md$Nr2f2[md$shf_status=="SHF_low"]; nhi <- md$Nr2f2[md$shf_status=="SHF_high"]
cat(sprintf("Nr2f2 gene   | %%pos low %.2f vs high %.2f | cliffs %.3f | p %.3g\n",
            mean(nlo>0), mean(nhi>0), cliffs(nlo,nhi), wilcox.test(nlo,nhi)$p.value))

# context: is SHF-low arterial venous-INTERMEDIATE between SHF-high arterial and the venous bed?
FetchData(sub, c("bed","shf_status","Venous_score","Nr2f2")) |>
  mutate(grp = case_when(
    bed=="Arterial_PAA" & shf_status=="SHF_high" ~ "Art_SHFhigh",
    bed=="Arterial_PAA" & shf_status=="SHF_low"  ~ "Art_SHFlow",
    bed=="Venous"                                ~ "Venous_bed")) |>
  filter(!is.na(grp)) |>
  group_by(grp) |> summarise(mean_venous = mean(Venous_score),
                             pct_Nr2f2 = mean(Nr2f2 > 0), n = n())
