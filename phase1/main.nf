#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * Phase 1 — Per-sample QC, doublet removal, normalization and clustering
 * =============================================================================
 * Runs once per sample, in parallel (capped per profile). Reads the sample
 * roster from params/samples.yml and processes each sample through:
 *
 *   CREATE_SEURAT_OBJECT -> VISUALIZE_QC -> FILTER_CELLS -> DETECT_DOUBLETS
 *      -> NORMALIZE_DATA -> DIM_REDUCTION_AND_CLUSTER
 *
 * DESIGN NOTES for THIS project (diverges deliberately from the mesp1 guide):
 *   * DOUBLET DETECTION added at step 05 — AFTER basic QC filtering (so the
 *     doublet simulation isn't polluted by empty/dying cells) and BEFORE
 *     normalization/integration (so doublets never enter SCTransform's noise
 *     model or get blended into real clusters where they can't be found).
 *     It is per-sample because doublets form within a single 10x capture;
 *     cross-sample "doublets" are meaningless.
 *   * CLUSTREE REMOVED. Phase 1 clusters are a QC diagnostic, not a deliverable
 *     (the analytical clustering happens on the integrated atlas in Phase 3 and
 *     in the caudal-arch subset in Phase 4a). Resolution-stability rigor is
 *     spent there, not on per-sample clusters we discard.
 *
 * INPUT CONTRACT: each sample's raw 10x triplet must exist at
 *   <raw_data_dir>/<id>/{barcodes,features,matrix}.* — produced by the
 *   standalone scripts/download_geo.sh. Run that FIRST:
 *      bash scripts/download_geo.sh
 *      nextflow run phase1/main.nf -profile conda,apple_silicon
 *
 * KEY OUTPUTS per sample (consumed by later phases):
 *   results/phase1/<id>/objects/05_seurat_singlets.rds   -> Phase 2a integration
 *   results/phase1/<id>/objects/07_seurat_clustered.rds  -> Phase 2 annotation
 *
 * SYNTAX NOTE: under Nextflow >=25 strict syntax, only DECLARATIONS (processes,
 * workflows, functions, includes) are allowed at the top level — all
 * statements (channel construction, variable defs) must live INSIDE workflow{}.
 * =============================================================================
 */

// ============================ PROCESSES ====================================

// 02 — Build a Seurat object from the 10x sparse matrix
process CREATE_SEURAT_OBJECT {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy'
  input:
    tuple val(meta), path(raw_dir)
  output:
    tuple val(meta), path("02_seurat_unfiltered.rds"), emit: obj
  script:
  """
  Rscript ${projectDir}/scripts/02_create_seurat_obj.R \\
    --raw_dir      ${raw_dir} \\
    --sample_id    ${meta.id} \\
    --genotype     ${meta.genotype} \\
    --timepoint    ${meta.timepoint} \\
    --series_id    ${meta.series_id} \\
    --mt_pattern   '${params.mt_pattern}' \\
    --out          02_seurat_unfiltered.rds
  """
}

// 03 — Compute QC metrics; emit violin/scatter plots (no filtering yet)
process VISUALIZE_QC {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy', pattern: "*.rds"
  publishDir { "${params.results_dir}/${meta.id}/plots" },   mode: 'copy', pattern: "*.pdf"
  input:
    tuple val(meta), path(obj)
  output:
    tuple val(meta), path("03_seurat_with_qc.rds"), emit: obj
    path "03_qc_plots.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/03_qc_visualize.R \\
    --in_rds         ${obj} \\
    --sample_id      ${meta.id} \\
    --min_features   ${meta.min_features} \\
    --max_features   ${meta.max_features} \\
    --max_percent_mt ${meta.max_percent_mt} \\
    --out_rds        03_seurat_with_qc.rds \\
    --out_pdf        03_qc_plots.pdf
  """
}

// 04 — Remove low-quality cells (nFeature + mt%); thresholds are per-sample.
//      This is the "first-pass QC" scDblFinder expects upstream of it.
process FILTER_CELLS {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy'
  input:
    tuple val(meta), path(obj)
  output:
    tuple val(meta), path("04_seurat_filtered.rds"), emit: obj
  script:
  """
  Rscript ${projectDir}/scripts/04_filter_cells.R \\
    --in_rds         ${obj} \\
    --min_features   ${meta.min_features} \\
    --max_features   ${meta.max_features} \\
    --max_percent_mt ${meta.max_percent_mt} \\
    --out            04_seurat_filtered.rds
  """
}

// 05 — Per-sample doublet detection (scDblFinder). See 05_detect_doublets.R for
//      the trajectory-safety rationale (random doublet generation by default).
process DETECT_DOUBLETS {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy', pattern: "*.rds"
  publishDir { "${params.results_dir}/${meta.id}/plots" },   mode: 'copy', pattern: "*.pdf"
  publishDir { "${params.results_dir}/${meta.id}/tables" },  mode: 'copy', pattern: "*.csv"
  input:
    tuple val(meta), path(obj)
  output:
    tuple val(meta), path("05_seurat_singlets.rds"), emit: obj   // -> Phase 2a
    path "05_doublet_summary.csv"
    path "05_doublet_plots.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/05_detect_doublets.R \\
    --in_rds      ${obj} \\
    --sample_id   ${meta.id} \\
    --seed        ${params.seed} \\
    --out_rds     05_seurat_singlets.rds \\
    --out_summary 05_doublet_summary.csv \\
    --out_pdf     05_doublet_plots.pdf
  """
}

// 06 — SCTransform normalization (regress percent.mt) on clean singlets
process NORMALIZE_DATA {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy'
  input:
    tuple val(meta), path(obj)
  output:
    tuple val(meta), path("06_seurat_normalized.rds"), emit: obj
  script:
  """
  Rscript ${projectDir}/scripts/06_normalize_data.R \\
    --in_rds ${obj} \\
    --n_hvgs ${params.n_hvgs} \\
    --seed   ${params.seed} \\
    --out    06_seurat_normalized.rds
  """
}

// 07 — PCA -> UMAP -> Louvain clustering (QC diagnostic; not the final clusters)
process DIM_REDUCTION_AND_CLUSTER {
  tag "${meta.id}"
  publishDir { "${params.results_dir}/${meta.id}/objects" }, mode: 'copy', pattern: "*.rds"
  publishDir { "${params.results_dir}/${meta.id}/plots" },   mode: 'copy', pattern: "*.pdf"
  input:
    tuple val(meta), path(obj)
  output:
    tuple val(meta), path("07_seurat_clustered.rds"), emit: obj   // -> Phase 2
    path "07_elbow.pdf"
    path "07_umap.pdf"
  script:
  """
  Rscript ${projectDir}/scripts/07_run_dim_reduction.R \\
    --in_rds      ${obj} \\
    --n_pcs       ${params.n_pcs} \\
    --cluster_res ${params.cluster_res} \\
    --seed        ${params.seed} \\
    --out_rds     07_seurat_clustered.rds \\
    --out_elbow   07_elbow.pdf \\
    --out_umap    07_umap.pdf
  """
}

// ============================ WORKFLOW =====================================
// All channel construction lives INSIDE the workflow block (strict-syntax rule:
// top-level statements are not allowed; only declarations).
workflow {

  // Read the sample roster (SnakeYAML ships with Nextflow; fully-qualified
  // class name is used because top-level `import` is also disallowed).
  def cfg       = new org.yaml.snakeyaml.Yaml().load(file(params.input_samples).text)
  def overrides = cfg.qc_overrides ?: [:]

  // Build the per-sample channel: tuple(meta, raw_dir).
  // Per-sample QC overrides win over the global defaults from nextflow.config.
  samples_ch = Channel
    .fromList(cfg.samples)
    .map { s ->
        def ov   = (overrides[s.id] ?: [:])
        def meta = [
          id            : s.id,
          genotype      : s.genotype,
          timepoint     : s.timepoint,
          series_id     : s.series_id,
          min_features  : (ov.min_features   != null ? ov.min_features   : params.min_features),
          max_features  : (ov.max_features   != null ? ov.max_features   : params.max_features),
          max_percent_mt: (ov.max_percent_mt != null ? ov.max_percent_mt : params.max_percent_mt)
        ]
        def raw_dir = file("${params.raw_data_dir}/${s.id}")
        tuple(meta, raw_dir)
    }

  // Wire the per-sample DAG.
  CREATE_SEURAT_OBJECT(samples_ch)
  VISUALIZE_QC(CREATE_SEURAT_OBJECT.out.obj)
  FILTER_CELLS(VISUALIZE_QC.out.obj)
  DETECT_DOUBLETS(FILTER_CELLS.out.obj)
  NORMALIZE_DATA(DETECT_DOUBLETS.out.obj)
  DIM_REDUCTION_AND_CLUSTER(NORMALIZE_DATA.out.obj)
}


