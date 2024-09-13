nextflow.enable.dsl=2

include { rmats_prep as rmats_prep_g1; rmats_prep as rmats_prep_g2; rmats_post } from "./rmats"

workflow {
  def required_params = ["bam_g1", "gtf", "readLength", "out_dir"]
  def missing_params = []

  required_params.each {
    param -> if (params[param] == null) missing_params.add(param)
  }

  if (missing_params.size() > 0) {
    exit 1, "Missing required parameters: ${missing_params.join(', ')}"
  }

  gtf = Channel
    .fromPath(params.gtf)
    .ifEmpty {exit 1, "No gtf file found in ${params.gtf}"}
    .first()

  bam_g1_files = Channel
    .fromPath(params.bam_g1)
    .ifEmpty {exit 1, "No bam files found in ${params.bam_g1}"}

  bam_g2_files = Channel.empty()
  if (params.bam_g2) {
    bam_g2_files = Channel
      .fromPath(params.bam_g2)
      .ifEmpty {exit 1, "No bam files found in ${params.bam_g2}"}
  }

  bam_g1_ids = bam_g1_files.reduce([]) { accum, v ->
    i = accum.size()
    accum.add("g1_${i}")
    accum
  }.flatten()
  bam_g2_ids = bam_g2_files.reduce([]) { accum, v ->
    i = accum.size()
    accum.add("g2_${i}")
    accum
  }.flatten()

  rmats_prep_g1(bam_g1_files, bam_g1_ids, gtf)
  rmats_prep_g2(bam_g2_files, bam_g2_ids, gtf)

  rmats_post(
    rmats_prep_g1.out.bam_name.toList(),
    rmats_prep_g2.out.bam_name.toList(),
    rmats_prep_g1.out.rmats.toList(),
    rmats_prep_g2.out.rmats.toList(),
    gtf
  )
}
