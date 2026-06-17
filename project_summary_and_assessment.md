# paa-dev-trajectory — Project Summary & Assessment

*WT scRNA-seq of microdissected cardiopharyngeal region, GSE158941. Four timepoints (E8.0, E8.25, E9.5, E10.5), Mesp1-lineage cells. Focus: SHF contribution to the 4th pharyngeal arch artery (PAA) endothelium.*

---

## 1. Biological objectives

The lab's in vivo work establishes that the second heart field (SHF) contributes >60% of 4th-arch-artery endothelium (Mef2c-AHF-Cre, Isl1-Cre tracing), and that conditional deletion of *Itga5*/*Fn1* (or *Tbx1*/*Vegfr3*) reduces SHF contribution — vessels still form but with a developmental **delay at E10.5**. The 4th PAA is uniquely sensitive. *Itga5/Itgb1/Fn1 are the perturbation handle, not identity markers.*

Three questions were brought to the scRNA-seq:

1. **Identify the non-SHF "compensatory" cells** that build the arch artery when SHF input drops.
2. **Find the genes driving the delay and the plexus → patent-artery remodeling.**
3. **Find the genes orchestrating incorporation of compensatory cells.**

---

## 2. Objectives scorecard

| Objective | Status | Basis |
|---|---|---|
| 1. Non-SHF compensatory cells | **Redirected — not resolvable in WT** | Three convergent null tests; no distinct population exists at baseline. Requires cKO-vs-WT. |
| 2. Remodeling / arterialization gene program | **In progress — method now sound, not yet executed** | Extractable from WT as a continuous arterialization trajectory; approach defined this session. |
| 3. Compensatory-cell incorporation genes | **Parked for cKO** | Logically downstream of Objective 1. |

The headline: **the WT atlas has largely done the job it can do.** It defined the SHF-derived arterial/PAA endothelium (the reference arm for any future perturbation experiment) and showed rigorously that the compensation question is a *perturbation* question, not a baseline one. Objective 2 is the remaining live WT deliverable.

---

## 3. Questions answered (with confidence)

| Question | Answer | Confidence |
|---|---|---|
| Is there a SHF-derived arterial/PAA endothelial population in the data? | **Yes** — the Arterial_PAA bed: Tmem100/Epas1/Gja5/Efnb2/Dll4/Hey1/Sox17, Cxcr4⁺ (pharyngeal). | High |
| Does its SHF fraction match the in vivo >60%? | **Yes** — ~70% SHF-score-positive (130/56 within the bed). | Medium-high — *score ≠ lineage* (see lessons). |
| Does the arterial/PAA endothelium express the integrin axis? | **Yes** — Itga5/Itgb1/Fn1 expressed there (reported as an outcome, never used to annotate). | High |
| Is there a distinct non-SHF compensatory *signature* in WT? | **No.** | High — convergent nulls. |
| Are the SHF-low arterial cells tip-like / posterior / venous in origin? | **No** (tip null; posterior a weak non-significant trend; venous null). | High |

---

## 4. Analytical arc — what was done, and its reproducibility status

### Captured by the Nextflow pipeline (reproducible, orchestrated)

| Step | Output | In Nextflow? |
|---|---|---|
| Per-sample QC → Seurat → scDblFinder doublet removal → SCTransform (diagnostic) | phase1 objects | **Yes** |
| LogNormalize + HVG → fastMNN integration | `02_integrated.rds` | **Yes** |
| Clustering, cell-cycle scoring, marker DE | `03`–`06` | **Yes** |
| Whole-embryo lineage annotation (16 lineages) | `07_annotated.rds` | **Yes** |
| Pijuan-Sala label transfer | `08_*` | **Yes** |

Everything through whole-embryo annotation is in the pipeline and reproducible end-to-end (phase1 ~5 min, phase2a ~25 min). Methodological decisions are locked (LogNormalize feeds fastMNN; SCTransform diagnostic-only; cell-cycle regression deliberately excluded for developmental data).

### Standalone R scripts — modular & individually reproducible, but NOT in Nextflow

These are the endothelial-focused analysis. Each uses the house style (optparse params, `set.seed(42)`, `out()` helper, documented caveats), so each is reproducible on its own — but they are run by hand, not orchestrated, and they are exploratory rather than a settled pipeline.

| Script | Purpose | Reproducible | Notes |
|---|---|---|---|
| `09_subset_endothelium.R` | Subset endothelium, re-HVG, integrate (`--integrate TRUE`) **or** uncorrected (`FALSE`), recluster | Yes | Produces two objects: `endo_only/` (integrated, 6 beds, cell-type) and `endo_traj/` (uncorrected, geometry for trajectory). **`FindClusters` is non-deterministic across runs — freeze the chosen realization.** |
| `10_shf_modules.R` | Freeze bed annotation; SHF/Tip/PostHox module scores; SHF-high/low split; hypothesis tests | Yes | Anti-double-dipping built in (test genes exclude the DE hits). |
| `11_trajectory.R` (v2) | Clean structural Slingshot trajectory | Yes | **Superseded** — see Lessons; the foreign-label issue makes its topology unreliable. |
| `11a_traj_native_markers.R` | Verify `endo_traj` composition; regenerate native-cluster markers | Yes | Confirmed 1049 endothelium-only cells; regenerated correct native markers. |

**Reproducibility weakness:** repeated stale-file / wrong-object incidents, because `09_*` filenames are identical across `endo_only/`, `endo_subset/`, and `endo_traj/`, and the scripts write bare relative filenames (designed for Nextflow `publishDir`, but run standalone they collide). This caused several detours (annotating the wrong cluster set, loading the wrong object interactively).

---

## 5. Key results

**Endothelial bed structure (integrated, 1049 cells, 6 beds):** Venous, Posterior_sinusoidal, RA_producing, **Arterial_PAA**, Angioblast, Venous_lymphatic. Dropping endocardium was essential — its EndoMT signal had been masking arterial structure.

**The Arterial_PAA bed is the SHF-derived PAA endothelium.** Highest SHF module score of all beds (~70% positive), arterial program complete, pharyngeal (Cxcr4⁺), integrin axis present. This is the central, solid positive result.

**No distinct compensatory population in WT.** Within the arterial bed, SHF-low cells (the candidate non-SHF/compensatory cells) were tested against three biologically-motivated hypotheses, all with independent (non-circular) gene sets:

- *Angiogenic tip / sprouting* — null (Cliff's δ 0.087, p 0.35; clean re-test still p 0.30).
- *Posterior-Hox positional identity* — weak, non-significant trend (δ 0.171, p 0.065, p_adj 0.13).
- *Venous origin (Nr2f2)* — null; venous score is actually *lowest* in SHF-low cells, Nr2f2 flat.

The convergence is itself the finding: the SHF-low arterial cells are the low tail of a continuous SHF-score distribution within transcriptionally uniform committed arterial endothelium — not a separate cell type. Compensation is a perturbation response with no standing WT population.

**Native uncorrected clusters (8, stage-structured)** — annotated on canonical identity only (integrin excluded): Angioblast_E80 (root), RA_producing_E80, Venous_E80, Posterior_sinusoidal_E80, Progenitor_E825, **Arterial_E825**, **Arterial_E95**, Venous_capillary_E105. The arterial/PAA lineage shows up as the stage series Arterial_E825 → Arterial_E95. *Flag:* the E10.5 cluster (Venous_capillary) is heavily ambient-hemoglobin-contaminated — drop Hb genes before any late-window work.

---

## 6. Methodological lessons (hard-won, worth keeping)

1. **Two embeddings, two jobs.** Integrated (fastMNN) for cell-type identity — stage mixed out, gives clean beds. Uncorrected (PCA) for trajectory — stage preserved, the developmental axis intact. Using one for both is the error.
2. **The uncorrected data cannot yield clean cell-type clusters at any resolution.** Stage is the dominant variance axis (real development + technical batch, confounded since each timepoint is its own sample). High resolution splits type × stage; low resolution merges to stage-islands; neither isolates cell type. Cell-type identity has exactly one clean source: the integrated beds.
3. **The data is a continuum.** Endothelial identity (arterial↔venous↔capillary) is a gradient; the only discrete structure is the stage sampling. This argues against cluster-centroid trajectory methods and for a continuous, score-ordered approach.
4. **Slingshot needs cluster labels coherent in the embedding it runs on.** Feeding integrated bed labels onto the uncorrected geometry produced a biologically false topology (arterial hanging off a venous trunk) — the labels' centroids landed in meaningless in-between positions. This, plus the stage gaps, is why the full-timecourse trajectory was unusable.
5. **SHF score ≠ SHF lineage.** The module score reads current program expression; SHF-low conflates true non-SHF cells with SHF-derived cells that downregulated the program. Only the lab's lineage tracing or a perturbation disambiguates — which is the structural reason WT can't close the compensation question.

---

## 7. Open leads

- **The arterialization / remodeling gene program (Objective 2) — the live WT analysis.** Build an integrin-free arterial *identity* score (Gja5/Tmem100/Sox17/Dll4/Hey1) and order cells along it continuously, either as the developmental series (Angioblast → Arterial_E825 → Arterial_E95) or within a cleaned late window. The genes moving along that axis are the arterialization program — the closest WT answer to "what drives plexus→artery remodeling." *Method now defined; not yet executed.*
- **Posterior-Hox weak trend.** δ 0.17, underpowered. Not a finding, but directionally consistent with arch-positional identity (the microdissection makes posterior Hox a *positional* signal, 4th/6th arch, not contamination). Hold lightly; revisit with more cells.
- **The cKO-vs-WT experiment.** The definitive route for Objectives 1 and 3. This WT atlas is the reference arm. Worth scoping now: genotypes (*Itga5*/*Fn1* cKO), timepoints (E9.5/E10.5 to capture the delay), cell numbers to power a rare-population contrast.

---

## 8. Overall assessment & recommended next moves

**What's solid:** the SHF-derived Arterial_PAA bed definition; the three convergent null tests (rigorous, independent, correctly interpreted); the methodological understanding of why this dataset behaves as it does.

**What's incomplete:** the remodeling gene program (path clear, not run); the trajectory (earlier versions compromised by foreign labels + stage gaps — the corrected continuous approach has not yet been executed).

**What's correctly closed:** the compensation question, as a WT question. The nulls are not a failure; they are a clean negative that redirects the project to the cKO.

**Recommended moves, in order:**

1. **Decide the scope of the WT endgame.** Either deliver Objective 2 (the continuous arterialization program) as the final WT result, or declare the WT atlas complete-as-reference and pivot. Given the method is now sound, executing Objective 2 is a contained, high-value finish.
2. **Harden reproducibility before any more analysis.** Fold `09`–`11a` into the Nextflow pipeline (or a documented sub-workflow) with explicit, per-run output directories, and freeze each subset object under a unique, collision-proof name. Most of this session's detours were file/object-identity confusion, not science.
3. **Fix the ambient-Hb contamination** in the late/E10.5 endothelial cells (exclude Hb genes from HVGs) before any late-window analysis.
4. **Scope the cKO-vs-WT design** as the experiment the WT nulls point to — this is where Objectives 1 and 3 actually get answered.

**One-line verdict:** the WT analysis succeeded at what WT can do — it defined the SHF-derived PAA endothelium and proved the compensation question needs perturbation data. The remaining WT win is the arterialization program; the bigger scientific payoff is the cKO this work now clearly motivates.