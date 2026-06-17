#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * Phase 2a — WT Mesp1-lineage integration → clustering → annotation
 * =============================================================================
 * Fan-in pipeline: the four phase1 singlet objects are merged, integrated with
 * fastMNN, clustered at a locked resolution, cell-cycle scored, marker-tested,
 * and annotated against a canonical lineage panel.
 *
 *   NORMALIZE -> INTEGRATE -> CLUSTER -> CELLCYCLE -> MARKERS -> ANNOTATE
 *   (ANNOTATE also consumes 05_markers_top.csv from MARKERS)
 *
 * INPUT CONTRACT: each sample's phase1 output must exist at
 *   results/phase1/<id>/objects/05_seurat_singlets.rds
 * The samplesheet 01_normalize.R consumes is GENERATED from params/samples.yml
 * at launch (same single source of truth as phase1 / download) — no hand-edited
 * paths. Run phase1 first, then:
 *   nextflow run phase2a/main.nf -profile conda,apple_silicon
 *
 * SYNTAX NOTE: under Nextflow >=25 strict syntax, only DECLARATIONS live at the
 * top level — all statements (roster parse, samplesheet build) go inside workflow{}.
 * =============================================================================
 */

// ============================ PROCESSES ====================================

// 01 — Merge the four singlet objects, log-normalise, select HVGs
process NORMALIZE {
  tag "normalize"
  publishDir "${params.outdir}/objects", mode: 'copy', pattern: "*.rds"
  input:
    path samplesheet
  output:
    path "01_merged.rds", emit: obj
  script:
  """
  Rscript ${projectDir}/scripts/01_normalize.R \\
    --samplesheet ${samplesheet} \\
    --n_hvgs      ${params.n_hvgs}
  """
}

// 02 — fastMNN integration (temporal merge order) + mixing-matrix diagnostic
process INTEGRATE {
  tag "integrate"
  publishDir "${params.outdir}/objects", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/plots",   mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"
  input:
    path obj
  output:
    path "02_integrated.rds", emit: obj
    path "02_mixing_baseline.csv"
    path "02_mixing_mnn.csv"
    path "02_mixing_comparison.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/02_integrate.R \\
    --input  ${obj} \\
    --n_dims ${params.n_dims}
  """
}

// 03 — kNN graph on MNN embedding, cluster at locked resolution
process CLUSTER {
  tag "cluster"
  publishDir "${params.outdir}/objects", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/plots",   mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"
  input:
    path obj
  output:
    path "03_clustered.rds", emit: obj
    path "03_umap_clusters.pdf"
    path "03_cluster_summary.csv"
  script:
  """
  Rscript ${projectDir}/scripts/03_cluster.R \\
    --input      ${obj} \\
    --n_dims     ${params.n_dims} \\
    --resolution ${params.resolution} \\
    --algorithm  ${params.algorithm}
  """
}

// 04 — Cell-cycle scoring (Tirosh, diagnostic only; no regression — see script)
process CELLCYCLE {
  tag "cellcycle"
  publishDir "${params.outdir}/objects", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/plots",   mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"
  input:
    path obj
  output:
    path "04_cc_scored.rds", emit: obj
    path "04_umap_phase.pdf"
    path "04_phase_composition.pdf"
    path "04_phase_by_cluster.csv"
  script:
  """
  Rscript ${projectDir}/scripts/04_cellcycle.R \\
    --input      ${obj} \\
    --resolution ${params.resolution}
  """
}

// 05 — Unbiased cluster DE (Wilcoxon, positive markers only)
process MARKERS {
  tag "markers"
  publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.csv"
  publishDir "${params.outdir}/plots",  mode: 'copy', pattern: "*.pdf"
  input:
    path obj
  output:
    path "05_markers_all.csv"
    path "05_markers_top.csv", emit: top
    path "05_dotplot_markers.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/05_markers.R \\
    --input      ${obj} \\
    --resolution ${params.resolution} \\
    --top_n      ${params.top_n_markers} \\
    --min_pct    ${params.min_pct} \\
    --logfc      ${params.logfc}
  """
}

// 06 — Canonical lineage panel: dotplot + per-lineage module scores
process ANNOTATE {
  tag "annotate"
  publishDir "${params.outdir}/plots",  mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/tables", mode: 'copy', pattern: "*.csv"
  input:
    path obj
    path top_markers
  output:
    path "06_dotplot_canonical.pdf"
    path "06_score_heatmap.pdf"
    path "06_module_scores.csv"
  script:
  """
  Rscript ${projectDir}/scripts/06_annotate.R \\
    --input      ${obj} \\
    --markers    ${top_markers} \\
    --resolution ${params.resolution}
  """
}

process FREEZE {
  tag "freeze"
  publishDir "${params.outdir}/objects", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/plots",   mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/tables",  mode: 'copy', pattern: "*.csv"
  input:
    path obj
  output:
    path "07_annotated.rds", emit: obj
    path "07_lineage_map.csv"
    path "07_umap_lineage.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/07_annotate_freeze.R --input ${obj} --resolution ${params.resolution}
  """
}

// ============================ WORKFLOW =====================================
workflow {

  // ── Roster = single source of truth (same SnakeYAML idiom as phase1) ───────
  def cfg = new org.yaml.snakeyaml.Yaml().load(file(params.input_samples).text)

  // ── Generate the samplesheet from the roster ───────────────────────────────
  // One row per sample, pointing at its phase1 singlet object. No hand-edited
  // paths: ids come from the roster, the rest is fixed by the phase1 layout.
  // Fail fast (naming the offending id) if a singlet is missing.
  def rows = cfg.samples.collect { s ->
      def rds = file("${params.singlets_dir}/${s.id}/objects/05_seurat_singlets.rds")
      if( !rds.exists() )
          error "Missing singlets for '${s.id}': ${rds}\n  -> run phase1 first."
      "${s.id},${rds},${s.timepoint}"
  }

  def samplesheet = file("${params.outdir}/samplesheet.csv")
  samplesheet.parent.mkdirs()
  samplesheet.text = (['sample_id,path,timepoint'] + rows).join('\n') + '\n'

  // ── Linear fan-in DAG ──────────────────────────────────────────────────────
  NORMALIZE( samplesheet )
  INTEGRATE( NORMALIZE.out.obj )
  CLUSTER(   INTEGRATE.out.obj )
  CELLCYCLE( CLUSTER.out.obj )
  MARKERS(   CELLCYCLE.out.obj )
  ANNOTATE(  CELLCYCLE.out.obj, MARKERS.out.top )
  FREEZE( CELLCYCLE.out.obj )
}
