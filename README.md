# Mesodermal-Derived Endothelial Clusters

Single-cell RNA-seq atlas of Mesp1-lineage cardiovascular progenitors resolving the endothelial compartment of the developing cardiopharyngeal field, with focus on the second heart field (SHF) contribution to fourth pharyngeal arch artery (PAA) endothelium.

---

## Biological Context

The fourth pharyngeal arch artery is the embryonic precursor of the aortic arch, a structure whose developmental failure underlies the cardiovascular malformations seen in 22q11.2 deletion syndrome (DiGeorge). In vivo lineage tracing (Mef2c-AHF-Cre, Isl1-Cre) shows that second heart field progenitors contribute >75% of the endothelial cells lining the fourth PAA, yet the molecular identity of those SHF-derived endothelial cells, what distinguishes them from non-SHF endothelium, and which genes orchestrate the plexus-to-artery remodeling transition remain unresolved.

This project constructs a whole-embryo single-cell atlas of WT Mesp1-lineage cells spanning four developmental timepoints (E8.0 → E10.5), subsets and re-resolves the endothelial compartment into named arterial/venous/progenitor beds, characterises the SHF contribution to the arterial bed, and defines the arterialization gene program active in the late developmental window.

**Three analytical objectives:**

| # | Objective | Status |
|---|---|---|
| 1 | Define SHF-derived PAA arterial endothelial cells; identify compensatory non-SHF population | **Complete — key result confirmed** |
| 2 | Genes driving the plexus→artery remodeling/arterialization transition | **In progress (script 04)** |
| 3 | Genes orchestrating compensatory-cell incorporation (WT → cKO comparison) | **Blocked — deferred to cKO-vs-WT project** |

---

## Data

| Field | Value |
|---|---|
| GEO accession | [GSE158941](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158941) |
| Organism | *Mus musculus* |
| Tissue | Microdissected cardiopharyngeal region |
| Lineage | Mesp1-Cre; GFP+ sorted cells |
| Timepoints | E8.0, E8.25, E9.5, E10.5 |
| Cells (post-QC) | ~22,000 whole-embryo; 1,049 endothelium-only |
| Platform | Illumina (10x Chromium) |
| Genotype | Wild-type |

---

## Repository Structure

```
mesodermal-derived-endothelial-clusters/
├── phase1/                     # Per-sample QC (Nextflow)
│   ├── main.nf
│   ├── nextflow.config
│   └── scripts/
│       ├── 02_create_seurat_obj.R
│       ├── 03_qc_visualize.R
│       ├── 04_filter_cells.R
│       ├── 05_detect_doublets.R
│       ├── 06_normalize_data.R
│       └── 07_run_dim_reduction.R
├── phase2/                     # Whole-embryo integration → annotation (Nextflow)
│   ├── main.nf
│   ├── nextflow.config
│   └── scripts/
│       ├── 01_normalize.R
│       ├── 02_integrate.R
│       ├── 03_cluster.R
│       ├── 04_cellcycle.R
│       ├── 05_markers.R
│       ├── 06_annotate.R
│       ├── 07_annotate_freeze.R
│       └── 08_label_transfer.R
├── phase3/                     # Per-lineage interrogation (Nextflow + standalone)
│   ├── main.nf                 # orchestrates 01 (both embeddings) + 02
│   ├── nextflow.config
│   ├── scripts/
│   │   ├── 01_subset_endothelium.R
│   │   ├── 02_shf_modules.R
│   │   ├── 03_trajectory.R
│   │   ├── 03a_traj_native_markers.R
│   │   └── 04_arterialization.R
│   └── exploratory/
│       └── Hypothesis_testing_SHF_v_nSHF.R
├── params/
│   └── samples.yml             # single source of truth for sample roster
├── scripts/
│   ├── download_geo.py
│   └── download_geo.sh
├── results/                    # gitignored — generated locally
├── raw_data/                   # gitignored
├── work/                       # gitignored (Nextflow scratch)
└── project_summary_and_assessment.md
```

---

## Pipeline Overview

### Phase 1 — Per-sample QC

One Seurat object per sample. QC metrics computed (nFeature, nCount, %MT), thresholds applied (500–7000 features, <5% MT), doublets detected (DoubletFinder), and per-sample normalization + dim reduction run for QC inspection.

```bash
nextflow run phase1/main.nf -profile conda,apple_silicon
```

### Phase 2 — Whole-Embryo Integration and Annotation

Four singlet objects merged, log-normalised, and integrated with fastMNN (temporal merge order E8.0→E8.25→E9.5→E10.5). Clustered at resolution 0.4 (locked from clustree sweep). Cell-cycle scored (Tirosh; no regression — documented in script 04). Unbiased markers + canonical lineage panel dotplot. Lineage labels frozen into a `lineage` metadata column (`07_annotate_freeze.R`). Reference label transfer from Pijuan-Sala E8.0–8.5 atlas resolves contested calls (script 08).

```bash
nextflow run phase2/main.nf -profile conda,apple_silicon
```

**Whole-embryo clusters (resolution 0.4, 16 lineages):**

| Cluster | Lineage | Key markers |
|---|---|---|
| 0 | Lateral_plate_mesoderm | Osr1, Foxf1, Fendrr |
| 1 | Second_heart_field | Sfrp5, Mab21l2 |
| 2 | Cranial_mesenchyme | Alx3 |
| 3 | Cranial_paraxial_mesoderm | Otx2, Zic2 |
| 4 | Branchiomeric_muscle_progenitor | Myf5, Tbx1, Isl1, Pitx2 |
| 5 | Somitic_mesoderm | Meox1, Pax3 |
| 6 | Allantois | Tbx4, posterior Hox |
| 7 | Caudal_mesoderm | Cdx1, Cdx4 |
| 8 | Endothelium | Pecam1, Cdh5, Etv2 |
| 9 | Primitive_erythroid | Hemgn, Hbb-y |
| 10 | Pharyngeal_mesenchyme | Pax9, Sostdc1 |
| 11 | Hematopoietic | Rac2, Fermt3 |
| 12 | Cardiomyocyte | Nppa, Myh6 |
| 13 | Endocardium | Klf2, Bmx |
| 14 | Mesenchyme | reference transfer 0.995 |
| 15 | Epicardium | Upk3b, Aldh1a1 |

### Phase 3 — Endothelial Compartment Interrogation

Scripts 01 and 02 are orchestrated by Nextflow. Scripts 03a and 04 run standalone.

```bash
# scripts 01 (both embeddings) + 02
nextflow run phase3/main.nf -profile conda,apple_silicon

# native cluster verification (standalone)
Rscript phase3/scripts/03a_traj_native_markers.R

# arterialization gene program (standalone)
Rscript phase3/scripts/04_arterialization.R
```

**Two parallel embeddings from script 01:**

| Embedding | Integration | Purpose | Output subdir |
|---|---|---|---|
| Integrated | fastMNN (timepoints mixed out) | Cell-type beds | `results/phase3/endothelium/integrated/` |
| Uncorrected | PCA (timepoints preserved) | Developmental ordering | `results/phase3/endothelium/uncorrected/` |

**Endothelial beds (fastMNN-integrated, resolution 0.6, 1,049 cells):**

| Subcluster | Bed | Key markers |
|---|---|---|
| 0 | Venous | Ackr1 |
| 1 | Posterior_sinusoidal | Hoxb5, Oit3 |
| 2 | RA_producing | Aldh1a2, Bmp2 |
| 3 | Arterial_PAA | Tmem100, Gja5, Sox17, Dll4, Hey1, Efnb2 |
| 4 | Angioblast | Etv2, Sp5 |
| 5 | Venous_lymphatic | Nr2f2, Lyve1, Calcrl |

---

## Key Results

### Objective 1 — SHF contribution to Arterial_PAA endothelium

Bed 3 (Arterial_PAA) carries the highest SHF module score across all endothelial beds (~70% SHF-positive by AddModuleScore thresholding). This matches the >75% SHF contribution observed in vivo by lineage tracing and directly identifies the SHF-derived arterial endothelial cells at single-cell resolution.

Three independent compensatory-signature tests in SHF-low Arterial_PAA cells (angiogenic tip program, posterior Hox program, venous/Nr2f2 program) were all null — no distinct compensatory non-SHF population exists in WT. The compensation question is correctly deferred to a cKO-vs-WT comparison (separate project).

### Trajectory analysis — negative result, empirically documented

Slingshot-based structural trajectory inference was tested and formally retired. Even a minimal 3-node restricted fit (Angioblast_E80 → Arterial_E825 → Arterial_E95, root anchored) inverted developmental order, producing Angioblast → Arterial_E95 → Arterial_E825. Root cause: developmental stage is the dominant axis of variance in the uncorrected embedding, and Euclidean distance between cluster centroids does not reliably track lineage adjacency even after restricting to a biologically unambiguous 3-node series. This is not a Slingshot-specific limitation — it applies to any cluster-chaining trajectory method on this embedding.

See `results/phase3/endothelium/uncorrected/03_lineage_structure.txt` and script `03_trajectory.R` header for full documentation.

### Objective 2 — Arterialization gene program (in progress)

Script `04_arterialization.R` defines a continuous arterialization axis using an integrin-free arterial identity score (Gja5, Tmem100, Sox17, Dll4, Hey1, Efnb2, Cxcr4, Epas1) across the late developmental window (E9.5/E10.5). Preliminary top program genes: Gja4 (ρ = 0.54), Nr2f2 (ρ = −0.40, correctly anti-correlated), Mecom (ρ = 0.38), Bmx (ρ = 0.35) — independent recovery of known arterial identity markers validates the scoring axis.

---

## Requirements

**R packages (conda env: paa-dev):**
Seurat ≥5, batchelor, BiocParallel, slingshot, tidyverse, glue, optparse, viridisLite, clustree, future, DoubletFinder

**Nextflow:** ≥24.04.0

**Hardware:** Developed on macOS M3 Max 64 GB. Phase 3 runs two 32 GB processes in parallel; set `process.maxForks = 1` in `phase3/nextflow.config` if memory is limited.

---

## Reproducing the Analysis

```bash
# 1. Clone
git clone git@github.com:drgideonobeng/mesodermal-derived-endothelial-clusters.git
cd mesodermal-derived-endothelial-clusters

# 2. Download raw data (GEO)
bash scripts/download_geo.sh   # downloads GSE158941 to raw_data/

# 3. Per-sample QC
nextflow run phase1/main.nf -profile conda,apple_silicon

# 4. Whole-embryo integration → annotation
nextflow run phase2/main.nf -profile conda,apple_silicon

# 5. Endothelial interrogation
nextflow run phase3/main.nf -profile conda,apple_silicon

# 6. Verification + arterialization (standalone)
conda run -n paa-dev Rscript phase3/scripts/03a_traj_native_markers.R
conda run -n paa-dev Rscript phase3/scripts/04_arterialization.R
```

All scripts are seeded (`set.seed(42)`) and parametrized via optparse. Run any script with `--help` for full option listing.

---

## Citation

Nomaru H, Liu Y, Kim HJ, et al. Single cell multi-omics analysis reveals the gene regulatory network underlying differentiation of the cardiopharyngeal mesoderm. *Nature Communications* (2021). PMID: 34789765

---

## Author

Gideon Obeng, PhD
Postdoctoral Associate | American Heart Association Fellow 
Rutgers New Jersey Medical School 
[github.com/drgideonobeng](https://github.com/drgideonobeng)
