#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * Phase 3 — per-lineage interrogation: endothelial compartment (scripts 01–02)
 * =============================================================================
 * Forks from the whole-embryo phase2 annotated object into two parallel
 * endothelial embeddings, then chains the integrated embedding into bed
 * annotation and SHF module scoring:
 *
 *   07_annotated.rds ─┬─ SUBSET_INTEGRATED  (01, fastMNN) ──→ SHF_MODULES (02)
 *                     └─ SUBSET_UNCORRECTED (01, PCA)      ──→ [03/04, separate]
 *
 * Design rationale (see script headers for full detail):
 *   SUBSET_INTEGRATED  — stage mixes out → clean cell-TYPE beds for identity
 *   SUBSET_UNCORRECTED — stage preserved → developmental axis for trajectory
 *   SHF_MODULES        — freezes the 6 bed labels; scores SHF/Tip/PostHox
 *                        programs; tests compensatory leads in Arterial_PAA
 *
 * Run:
 *   nextflow run phase3/main.nf -profile conda,apple_silicon
 * =============================================================================
 */

// ============================ PROCESSES =====================================

// 01a — Endothelial subset, fastMNN-integrated (bed clustering substrate)
process SUBSET_INTEGRATED {
  tag "subset_integrated"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.csv"
  input:
    path annotated_obj
  output:
    path "01_endo_subset.rds",           emit: obj
    path "01_endo_umap_overlays.pdf"
    path "01_endo_bed_dotplot.pdf"
    path "01_endo_subcluster_markers.csv"
    path "01_endo_subcluster_top.csv"
    path "01_endo_subcluster_summary.csv"
  script:
  """
  Rscript ${projectDir}/scripts/01_subset_endothelium.R \\
    --input           ${annotated_obj} \\
    --lineages        ${params.lineages_subset} \\
    --n_hvgs          ${params.n_hvgs} \\
    --n_dims          ${params.n_dims} \\
    --resolution      ${params.resolution} \\
    --k               ${params.k_mnn} \\
    --min_cells_batch ${params.min_cells_batch} \\
    --integrate       TRUE \\
    --lineage         ${params.lineage} \\
    --outdir          .
  """
}

// 01b — Endothelial subset, uncorrected PCA (trajectory/ordering geometry)
// Produces a separate embedding of the same cells; output staged to uncorrected/
// for downstream 03_trajectory and 04_arterialization (not wired in this module).
process SUBSET_UNCORRECTED {
  tag "subset_uncorrected"
  publishDir "${params.outdir}/${params.lineage}/uncorrected", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/${params.lineage}/uncorrected", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/${params.lineage}/uncorrected", mode: 'copy', pattern: "*.csv"
  input:
    path annotated_obj
  output:
    path "01_endo_subset.rds",           emit: obj
    path "01_endo_umap_overlays.pdf"
    path "01_endo_bed_dotplot.pdf"
    path "01_endo_subcluster_markers.csv"
    path "01_endo_subcluster_top.csv"
    path "01_endo_subcluster_summary.csv"
  script:
  """
  Rscript ${projectDir}/scripts/01_subset_endothelium.R \\
    --input           ${annotated_obj} \\
    --lineages        ${params.lineages_subset} \\
    --n_hvgs          ${params.n_hvgs} \\
    --n_dims          ${params.n_dims} \\
    --resolution      ${params.resolution} \\
    --k               ${params.k_mnn} \\
    --min_cells_batch ${params.min_cells_batch} \\
    --integrate       FALSE \\
    --lineage         ${params.lineage} \\
    --outdir          .
  """
}

// 02 — Freeze bed annotation; score SHF/Tip/PostHox programs;
//      test compensatory leads in the Arterial_PAA SHF split.
//      Consumes the INTEGRATED object only (beds are not defined in uncorrected).
process SHF_MODULES {
  tag "shf_modules"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.rds"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/${params.lineage}/integrated", mode: 'copy', pattern: "*.csv"
  input:
    path integrated_obj
  output:
    path "02_endo_annotated.rds", emit: obj
    path "02_module_violins.pdf"
    path "02_module_tests.csv"
  script:
  """
  Rscript ${projectDir}/scripts/02_shf_modules.R \\
    --input     ${integrated_obj} \\
    --focus_bed ${params.focus_bed} \\
    --lineage   ${params.lineage} \\
    --outdir    .
  """
}

// ============================ WORKFLOW =====================================
workflow {

  // Phase2 annotated object — single upstream input, no samplesheet needed.
  // Guard: fail fast rather than producing an empty or corrupt downstream object.
  def annotated_obj = file(params.annotated_obj)
  if ( !annotated_obj.exists() )
      error "Phase2 annotated object not found: ${annotated_obj}\n  -> run phase2 first."

  // Fork: both embeddings run in parallel from the same input
  SUBSET_INTEGRATED(  annotated_obj )
  SUBSET_UNCORRECTED( annotated_obj )

  // Chain: SHF scoring waits for the integrated embedding only
  SHF_MODULES( SUBSET_INTEGRATED.out.obj )
}
