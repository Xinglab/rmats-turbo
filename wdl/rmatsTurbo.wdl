## version 1.0
## This workflow "rMATS_turbo" takes BAM files as input.

workflow rMATS_turbo {
  Array[File] bam_g1
  Array[File] bam_g2
  File gtf
  String pairedORsingle = "paired"
  Int readLength
  Int nthread = 1
  String out_dir
  String lib_type = "fr-unstranded"
  Boolean variable_readLength = false
  Int? anchorLength
  Int tstat = 1
  Float cstat = 0.0001
  Boolean statoff = false
  Boolean paired_stats = false
  Boolean novelSS = false
  Int mil = 50
  Int mel = 500
  Boolean allow_clipping = false

  Int machine_mem_gb = 4
  Int disk_space_gb = 20
  Boolean use_ssd = false
  String rmats_version = "v4.1.2"

  scatter (i in range(length(bam_g1))) {
    call rmats_pre as rmats_pre1 {
      input:
      bam = bam_g1[i],
      bam_id = "g1_"+i,
      gtf = gtf,
      pairedORsingle = pairedORsingle,
      readLength = readLength,
      out_dir = out_dir,
      lib_type = lib_type,
      variable_readLength = variable_readLength,
      anchorLength = anchorLength,
      novelSS = novelSS,
      mil = mil,
      mel = mel,
      machine_mem_gb = machine_mem_gb,
      disk_space_gb = disk_space_gb,
      use_ssd = use_ssd,
      rmats_version = rmats_version,
      allow_clipping = allow_clipping
    }
  }

  scatter (i in range(length(bam_g2))) {
    call rmats_pre as rmats_pre2 {
      input:
      bam = bam_g2[i],
      bam_id = "g2_"+i,
      gtf = gtf,
      pairedORsingle = pairedORsingle,
      readLength = readLength,
      out_dir = out_dir,
      lib_type = lib_type,
      variable_readLength = variable_readLength,
      anchorLength = anchorLength,
      novelSS = novelSS,
      mil = mil,
      mel = mel,
      machine_mem_gb = machine_mem_gb,
      disk_space_gb = disk_space_gb,
      use_ssd = use_ssd,
      rmats_version = rmats_version,
      allow_clipping = allow_clipping
    }
  }

  call rmats_post {
    input:
    bam_name_g1 = rmats_pre1.bam_name,
    bam_name_g2 = rmats_pre2.bam_name,
    gtf = gtf,
    readLength = readLength,
    nthread = nthread,
    out_dir = out_dir,
    input_rmats1 = rmats_pre1.out_rmats,
    input_rmats2 = rmats_pre2.out_rmats,
    anchorLength = anchorLength,
    tstat = tstat,
    cstat = cstat,
    statoff = statoff,
    paired_stats = paired_stats,
    novelSS = novelSS,
    mil = mil,
    mel = mel,
    machine_mem_gb = machine_mem_gb,
    disk_space_gb = disk_space_gb,
    use_ssd = use_ssd,
    rmats_version = rmats_version
  }

  output {
    Array[File] read_outcomes = flatten([rmats_pre1.read_outcome, rmats_pre2.read_outcome])
    File out_tar = rmats_post.out_tar
  }

  meta {
    author: "Qian Liu"
    email : "Qian.Liu@RoswellPark.org"
    description: "WDL workflow on AnVIL for rMATS turbo v4.1.2 developed in Dr. Yi Xing's lab"
  }
}

task rmats_pre {
  File bam
  String bam_id
  File gtf
  String pairedORsingle
  Int readLength
  String out_dir
  String lib_type
  Boolean variable_readLength
  Int? anchorLength
  Boolean novelSS
  Int mil
  Int mel
  Boolean allow_clipping

  String variable_readLength_opt = if variable_readLength then "--variable-read-length" else ""
  String anchorLength_opt = if defined(anchorLength) then "--anchorLength" else ""
  String novelSS_opt = if novelSS then "--novelSS" else ""
  String mil_opt = if novelSS then "--mil" else ""
  String? mil_val = if novelSS then mil else ""
  String mel_opt = if novelSS then "--mel" else ""
  String? mel_val = if novelSS then mel else ""
  String allow_clipping_opt = if defined(allow_clipping) then "--allow-clipping" else ""

  Int machine_mem_gb
  Int disk_space_gb
  Boolean use_ssd
  String rmats_version

  command {
    echo ${bam} > prep.txt
    python /rmats/rmats.py --b1 prep.txt --gtf ${gtf} -t ${pairedORsingle} --readLength ${readLength} --nthread 1 --od ${out_dir} --tmp tmp_output_prep_${bam_id} --task prep --libType ${lib_type} ${variable_readLength_opt} ${anchorLength_opt} ${anchorLength} ${novelSS_opt} ${mil_opt} ${mil_val} ${mel_opt} ${mel_val} ${allow_clipping_opt}
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
  Float cstat
  Boolean statoff
  Boolean paired_stats
  Boolean novelSS
  Int mil
  Int mel

  String anchorLength_opt = if defined(anchorLength) then "--anchorLength" else ""
  String cstat_opt = if paired_stats then "" else "--cstat"
  String cstat_val = if paired_stats then "" else cstat
  String statoff_opt = if statoff then "--statoff" else ""
  String paired_stats_opt = if paired_stats then "--paired-stats" else ""
  String novelSS_opt = if novelSS then "--novelSS" else ""
  String mil_opt = if novelSS then "--mil" else ""
  String? mil_val = if novelSS then mil else ""
  String mel_opt = if novelSS then "--mel" else ""
  String? mel_val = if novelSS then mel else ""

  Int machine_mem_gb
  Int disk_space_gb
  Boolean use_ssd
  String rmats_version

  Array[File] rmats1 = flatten(input_rmats1)
  Array[File] rmats2 = flatten(input_rmats2)
  Array[File] rmats = flatten([rmats1, rmats2])

  command {
    mkdir fd_rmats
    for file in ${sep=" " rmats}; do fn=`basename $file`; sed 's/.*\///g' $file > fd_rmats/$fn; done
    echo ${sep="," bam_name_g1} > bam_g1.txt
    echo ${sep="," bam_name_g2} > bam_g2.txt
    python /rmats/rmats.py --b1 bam_g1.txt --b2 bam_g2.txt --gtf ${gtf} --readLength ${readLength} --nthread ${nthread} --od ${out_dir} --tmp fd_rmats --task post ${anchorLength_opt} ${anchorLength} --tstat ${tstat} ${cstat_opt} ${cstat_val} ${statoff_opt} ${paired_stats_opt} ${novelSS_opt} ${mil_opt} ${mil_val} ${mel_opt} ${mel_val}
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

  parameter_meta {
    out_dir: "The directory for final output."
    type: "Type of read used in the analysis: either 'paired' for paired-end data or 'single' for single-end data. (Default: paired)"
    lib_type: "{fr-unstranded,fr-firststrand,fr-secondstrand} Library type. Use 'fr-firststrand' or 'fr-secondstrand' for strand-specific data. (Default: fr-unstranded)"
    readLength: "The length of each read."
    variable_readLength: "{true, false} Allow reads with lengths that differ from --readLength to be processed. --readLength will still    be used to determine IncFormLen and SkipFormLen. (Default: false)"
    anchorLength: "The anchor length. (Default: 1)"
    nthread: "The number of threads. The optimal number of threads should be equal to the number of CPU cores. (Default: 1)"
    tstat: "The number of threads for the statistical model. (Default: 1)"
    cstat: "The cutoff splicing difference. The cutoff used in the null hypothesis test for differential splicing. Valid: 0 <= cutoff < 1. Does not apply to the paired stats model. (Default: 0.0001 for 0.01% difference)"
    statoff: "{true, false} If skip the statistical analysis. (Default: false)"
    paired_stats: "{true, false} Use the paired stats model. (Default: false)"
    novelSS: "{true, false} Enable detection of novel splice sites (unannotated splice sites). (Default: false)"
    mil: "Minimum Intron Length. Only impacts --novelSS behavior. (Default: 50)"
    mel: "Maximum Exon Length. Only impacts --novelSS behavior. (Default: 500)"
  }
}
