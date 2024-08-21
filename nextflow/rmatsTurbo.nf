nextflow.enable.dsl=2

include { rmats_prep as rmats_prep_g1; rmats_prep as rmats_prep_g2; rmats_post } from "./rmats"

workflow {
  def required_params = ["bam_g1", "bam_g2", "gtf", "readLength", "out_dir"]
  def missing_params = []

  required_params.each {
    param -> if (params[param] == null) missing_params.add(param)
  }

  if (missing_params.size() > 0) {
    exit 1, "Missing required parameters: ${missing_params.join(', ')}"
  }

  bam_g1_ch = Channel
    .fromPath(params.bam_g1)
    .ifEmpty {exit 1, "No bam files found in ${params.bam_g1}"}
    .map {
      bam -> tuple(bam.simpleName, bam)
    }
  bam_g2_ch = Channel
    .fromPath(params.bam_g2)
    .ifEmpty {exit 1, "No bam files found in ${params.bam_g2}"}
    .map {
      bam -> tuple(bam.simpleName, bam)
    }

  gtf_ch = Channel
    .fromPath(params.gtf)
    .ifEmpty {exit 1, "No gtf file found in ${params.gtf}"}
  rmats_prep_g1(bam_g1_ch, gtf_ch, "g1")
  rmats_prep_g2(bam_g2_ch, gtf_ch, "g2")
  rmats_post(
    bam_g1_ch.map {name, bam -> bam}.collect(),
    bam_g2_ch.map {name, bam -> bam}.collect(),
    rmats_prep_g1.out.rmat.collect(),
    rmats_prep_g2.out.rmat.collect(),
    gtf_ch
  )
}
