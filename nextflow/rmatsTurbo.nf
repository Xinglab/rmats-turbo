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

  bam_g1_indices = bam_g1_files.reduce([]) { accum, v ->
    i = accum.size()
    accum.add(i)
    accum
  }.flatten()
  bam_g2_indices = bam_g2_files.reduce([]) { accum, v ->
    i = accum.size()
    accum.add(i)
    accum
  }.flatten()

  rmats_prep_g1(bam_g1_files, bam_g1_indices, gtf, "g1")
  rmats_prep_g2(bam_g2_files, bam_g2_indices, gtf, "g2")

  def get_prep_outputs = multiMapCriteria { sorted_runs ->
    rmats: {
      results = []
      for (run in sorted_runs) {
        results.add(run[1])
      }
      results
    }()
    bam_name: {
      results = []
      for (run in sorted_runs) {
        results.add(run[3])
      }
      results
    }()
  }

  // Sort prep output by bam_i
  sorted_prep_g1 = rmats_prep_g1.out
    .toSortedList{ a, b -> a[0] <=> b[0] }
    .multiMap(get_prep_outputs)

  sorted_prep_g2 = rmats_prep_g2.out
    .toSortedList{ a, b -> a[0] <=> b[0] }
    .multiMap(get_prep_outputs)

  rmats_post(
    sorted_prep_g1.bam_name,
    sorted_prep_g2.bam_name,
    sorted_prep_g1.rmats,
    sorted_prep_g2.rmats,
    gtf
  )
}
