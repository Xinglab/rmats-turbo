version 1.0

workflow rMATS_turbo {
  input {
    Array[File] bam_g1
    Array[File] bam_g2
    File gtf
    Boolean is_single_end = false
    Int readLength
    Int nthread = 1
    String out_dir
    String lib_type = "fr-unstranded"
    Boolean variable_read_length = false
    Int? anchorLength
    Int tstat = 1
    String cstat = "0.0001"
    Boolean statoff = false
    Boolean paired_stats = false
    Boolean darts_model = false
    String darts_cutoff = "0.05"
    Boolean novelSS = false
    Int mil = 50
    Int mel = 500
    Boolean allow_clipping = false
    Boolean individual_counts = false
    Int machine_mem_gb = 4
    Int disk_space_gb = 20
    Boolean use_ssd = false
    String rmats_version = "v4.3.0"
  }

  scatter (i in range(length(bam_g1))) {
    call rmats_prep as rmats_prep1 {
      input:
      bam = bam_g1[i],
      bam_id = "g1_"+i,
      gtf = gtf,
      is_single_end = is_single_end,
      readLength = readLength,
      out_dir = out_dir,
      lib_type = lib_type,
      variable_read_length = variable_read_length,
      anchorLength = anchorLength,
      novelSS = novelSS,
      mil = mil,
      mel = mel,
      allow_clipping = allow_clipping,
      machine_mem_gb = machine_mem_gb,
      disk_space_gb = disk_space_gb,
      use_ssd = use_ssd,
      rmats_version = rmats_version
    }
  }

  scatter (i in range(length(bam_g2))) {
    call rmats_prep as rmats_prep2 {
      input:
      bam = bam_g2[i],
      bam_id = "g2_"+i,
      gtf = gtf,
      is_single_end = is_single_end,
      readLength = readLength,
      out_dir = out_dir,
      lib_type = lib_type,
      variable_read_length = variable_read_length,
      anchorLength = anchorLength,
      novelSS = novelSS,
      mil = mil,
      mel = mel,
      allow_clipping = allow_clipping,
      machine_mem_gb = machine_mem_gb,
      disk_space_gb = disk_space_gb,
      use_ssd = use_ssd,
      rmats_version = rmats_version
    }
  }

  call rmats_post {
    input:
    bam_name_g1 = rmats_prep1.bam_name,
    bam_name_g2 = rmats_prep2.bam_name,
    gtf = gtf,
    readLength = readLength,
    nthread = nthread,
    out_dir = out_dir,
    input_rmats1 = rmats_prep1.out_rmats,
    input_rmats2 = rmats_prep2.out_rmats,
    anchorLength = anchorLength,
    tstat = tstat,
    cstat = cstat,
    statoff = statoff,
    paired_stats = paired_stats,
    darts_model = darts_model,
    darts_cutoff = darts_cutoff,
    novelSS = novelSS,
    mil = mil,
    mel = mel,
    individual_counts = individual_counts,
    machine_mem_gb = machine_mem_gb,
    disk_space_gb = disk_space_gb,
    use_ssd = use_ssd,
    rmats_version = rmats_version
  }

  output {
    Array[File] read_outcomes = flatten([rmats_prep1.read_outcome, rmats_prep2.read_outcome])
    File out_tar = rmats_post.out_tar
  }

  meta {
    author: "Qian Liu"
    email : "Qian.Liu@RoswellPark.org"
    description: "WDL workflow on AnVIL for rMATS turbo v4.3.0 developed in Dr. Yi Xing's lab"
  }
}

task rmats_prep {
  input {
    File bam
    String bam_id
    File gtf
    Boolean is_single_end
    Int readLength
    String out_dir
    String lib_type
    Boolean variable_read_length
    Int? anchorLength
    Boolean novelSS
    Int mil
    Int mel
    Boolean allow_clipping

    Int machine_mem_gb
    Int disk_space_gb
    Boolean use_ssd
    String rmats_version
  }

  String read_type_value = if is_single_end then "single" else "paired"
  String variable_read_length_opt = if variable_read_length then "--variable-read-length" else ""
  String anchorLength_opt = if defined(anchorLength) then "--anchorLength" else ""
  String novelSS_opt = if novelSS then "--novelSS" else ""
  String mil_opt = if novelSS then "--mil" else ""
  String mil_val = if novelSS then mil else ""
  String mel_opt = if novelSS then "--mel" else ""
  String mel_val = if novelSS then mel else ""
  String allow_clipping_opt = if allow_clipping then "--allow-clipping" else ""

  command {
    echo ${bam} > prep.txt
    python /rmats/rmats.py --b1 prep.txt --gtf ${gtf} -t ${read_type_value} --readLength ${readLength} --nthread 1 --od ${out_dir} --tmp tmp_output_prep_${bam_id} --task prep --libType ${lib_type} ${variable_read_length_opt} ${anchorLength_opt} ${anchorLength} ${novelSS_opt} ${mil_opt} ${mil_val} ${mel_opt} ${mel_val} ${allow_clipping_opt}
    mkdir outfd
    python /rmats/cp_with_prefix.py prep_${bam_id}_ outfd tmp_output_prep_${bam_id}/*.rmats
    cp tmp_output_prep_${bam_id}/*read_outcomes_by_bam.txt prep_${bam_id}_read_outcomes_by_bam.txt
  }

  runtime {
    docker: "xinglab/rmats:" + rmats_version
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
  }

  output {
    Array[File] out_rmats = glob("outfd/*.rmats")
    File read_outcome = "prep_${bam_id}_read_outcomes_by_bam.txt"
    String bam_name = basename("${bam}")
  }
}

task rmats_post {
  input {
    Array[String] bam_name_g1
    Array[String] bam_name_g2
    File gtf
    Int readLength
    Int nthread
    String out_dir
    Array[Array[File]] input_rmats1
    Array[Array[File]] input_rmats2
    Int? anchorLength
    Int tstat
    String cstat
    Boolean statoff
    Boolean paired_stats
    Boolean darts_model
    String darts_cutoff
    Boolean novelSS
    Int mil
    Int mel
    Boolean individual_counts

    Int machine_mem_gb
    Int disk_space_gb
    Boolean use_ssd
    String rmats_version
  }

  String anchorLength_opt = if defined(anchorLength) then "--anchorLength" else ""
  Boolean is_default_stats = (!paired_stats) && (!darts_model)
  String cstat_opt = if is_default_stats then "--cstat" else ""
  String cstat_val = if is_default_stats then cstat else ""
  String statoff_opt = if statoff then "--statoff" else ""
  String paired_stats_opt = if paired_stats then "--paired-stats" else ""
  String darts_model_opt = if darts_model then "--darts-model" else ""
  String darts_cutoff_opt = if darts_model then "--darts-cutoff" else ""
  String darts_cutoff_val = if darts_model then darts_cutoff else ""
  String novelSS_opt = if novelSS then "--novelSS" else ""
  String mil_opt = if novelSS then "--mil" else ""
  String mil_val = if novelSS then mil else ""
  String mel_opt = if novelSS then "--mel" else ""
  String mel_val = if novelSS then mel else ""
  String individual_counts_opt = if individual_counts then "--individual-counts" else ""

  Array[File] rmats1 = flatten(input_rmats1)
  Array[File] rmats2 = flatten(input_rmats2)
  Array[File] rmats = flatten([rmats1, rmats2])

  command {
    mkdir fd_rmats
    for file in ${sep=" " rmats}; do fn=`basename $file`; sed 's/.*\///g' $file > fd_rmats/$fn; done
    echo ${sep="," bam_name_g1} > bam_g1.txt
    echo ${sep="," bam_name_g2} > bam_g2.txt
    python /rmats/rmats.py --b1 bam_g1.txt --b2 bam_g2.txt --gtf ${gtf} --readLength ${readLength} --nthread ${nthread} --od ${out_dir} --tmp fd_rmats --task post ${anchorLength_opt} ${anchorLength} --tstat ${tstat} ${cstat_opt} ${cstat_val} ${statoff_opt} ${paired_stats_opt} ${darts_model_opt} ${darts_cutoff_opt} ${darts_cutoff_val} ${novelSS_opt} ${mil_opt} ${mil_val} ${mel_opt} ${mel_val} ${individual_counts_opt}
    tar czf ${out_dir}.tar.gz ${out_dir}
  }

  runtime {
    docker: "xinglab/rmats:" + rmats_version
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
    cpu: nthread
  }

  output {
    File out_tar = "${out_dir}.tar.gz"
  }
}
