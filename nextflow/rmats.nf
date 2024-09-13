process rmats_prep {
  tag "rmats_prep"
  label "mid_memory"
  publishDir path: "${params.publishDir}/read_outcomes", mode: 'copy', pattern: "prep_*_read_outcomes_by_bam.txt"

  input:
    path bam
    val bam_id
    path gtf

  output:
    path "outfd/*.rmats", emit: rmats
    path "prep_${bam_id}_read_outcomes_by_bam.txt"
    val bam_name, emit: bam_name

  script:
    bam_name = bam.getName()
    read_type_value = params.is_single_end ? "single" : "paired"
    variable_read_length_opt = params.variable_read_length ? "--variable-read-length" : ""
    anchorLength_opt = params.anchorLength ? "--anchorLength ${params.anchorLength}" : ""
    novelSS_opt = params.novelSS ? "--novelSS" : ""
    mil_opt = (params.novelSS && params.mil) ? "--mil ${params.mil}" : ""
    mel_opt = (params.novelSS && params.mel) ? "--mel ${params.mel}" : ""
    allow_clipping_opt = params.allow_clipping ? "--allow-clipping" : ""

  """
  echo ${bam} > prep.txt

  python /rmats/rmats.py \
    --b1 prep.txt \
    --gtf ${gtf} \
    -t ${read_type_value} \
    --readLength ${params.readLength} \
    --nthread 1 \
    --od ${params.out_dir} \
    --tmp tmp_output_prep_${bam_id} \
    --task prep \
    --libType ${params.lib_type} \
    ${variable_read_length_opt} \
    ${anchorLength_opt} \
    ${novelSS_opt} \
    ${mil_opt} ${mel_opt} \
    ${allow_clipping_opt}

    mkdir outfd

    python /rmats/cp_with_prefix.py prep_${bam_id}_ outfd tmp_output_prep_${bam_id}/*.rmats

    cp tmp_output_prep_${bam_id}/*read_outcomes_by_bam.txt prep_${bam_id}_read_outcomes_by_bam.txt
  """
}

process rmats_post {
  tag "rmats_post"
  label "tera_memory"
  publishDir "${params.publishDir}", mode: 'copy'

  input:
    val bams_g1
    val bams_g2
    path rmats_g1
    path rmats_g2
    path gtf
  output:
    path "${params.out_dir}.tar.gz"

  script:
  has_g2 = bams_g2.size() > 0
  b2_opt = has_g2 ? "--b2" : ""
  b2_val = has_g2 ? "bam_g2.txt" : ""

  anchorLength_opt = params.anchorLength ? "--anchorLength ${params.anchorLength}" : ""
  is_default_stats = (!params.paired_stats) && (!params.darts_model)
  cstat_opt = is_default_stats ? "--cstat ${params.cstat}" : ""
  statoff_opt = params.statoff ? "--statoff" : ""
  paired_stats_opt = params.paired_stats ? "--paired-stats" : ""
  darts_model_opt = params.darts_model ? "--darts-model" : ""
  darts_cutoff_opt = params.darts_model ? "--darts-cutoff ${params.darts_cutoff}" : ""
  novelSS_opt = params.novelSS ? "--novelSS" : ""
  mil_opt = (params.novelSS && params.mil) ? "--mil ${params.mil}" : ""
  mel_opt = (params.novelSS && params.mel) ? "--mel ${params.mel}" : ""
  individual_counts_opt = params.individual_counts ? "--individual-counts" : ""

  rmats1 = rmats_g1.flatten()
  rmats2 = rmats_g2.flatten()
  rmats = rmats1 + rmats2

  """
  mkdir fd_rmats

  for file in ${rmats.join(' ')}
  do
    fn=\$(basename \$file)
    sed 's/.*\\///g' \$file > fd_rmats/\$fn
  done

  echo ${bams_g1.join(',')} > bam_g1.txt
  echo ${bams_g2.join(',')} > bam_g2.txt

  python /rmats/rmats.py \
    --b1 bam_g1.txt \
    ${b2_opt} ${b2_val} \
    --gtf ${gtf} \
    --readLength ${params.readLength} \
    --nthread ${params.nthread} \
    --od ${params.out_dir} \
    --tmp fd_rmats \
    --task post \
    ${anchorLength_opt} \
    --tstat ${params.tstat} \
    ${cstat_opt} \
    ${statoff_opt} \
    ${paired_stats_opt} \
    ${darts_model_opt} \
    ${darts_cutoff_opt} \
    ${novelSS_opt} \
    ${mil_opt} \
    ${mel_opt} \
    ${individual_counts_opt}

    tar czf ${params.out_dir}.tar.gz ${params.out_dir}
  """
}
