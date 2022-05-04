rmats_docker <- "xinglab/rmats:v4.1.2-no-entrypoint"

get_array_type <- function(item_type) {
    return(list(type = "array", items = item_type))
}
get_2d_array_type <- function(item_type) {
    return(get_array_type(get_array_type(item_type)))
}

## Workflow inputs
wf_bam_g1_input <- InputParam(id = "wf_bam_g1", type = "File[]", position = -1)
wf_bam_g2_input <- InputParam(id = "wf_bam_g2", type = "File[]", position = -1)
wf_gtf_input <- InputParam(id = "wf_gtf", type = "File", position = -1)
wf_pairedORsingle_input <- InputParam(id = "wf_pairedORsingle", type = "string?", default = "paired",
                                      position = -1)
wf_readLength_input <- InputParam(id = "wf_readLength", type = "int", position = -1)
wf_nthread_input <- InputParam(id = "wf_nthread", type = "int?", default = 1L, position = -1)
wf_out_dir_input <- InputParam(id = "wf_out_dir", type = "string", position = -1)
wf_lib_type_input <- InputParam(id = "wf_lib_type", type = "string?", default = "fr-unstranded",
                                position = -1)
wf_variable_readLength_input <- InputParam(id = "wf_variable_readLength", type = "boolean?",
                                           default = FALSE, position = -1)
wf_anchorLength_input <- InputParam(id = "wf_anchorLength", type = "int?", position = -1)
wf_tstat_input <- InputParam(id = "wf_tstat", type = "int?", default = 1L, position = -1)
wf_cstat_input <- InputParam(id = "wf_cstat", type = "float?", default = 0.0001, position = -1)
wf_statoff_input <- InputParam(id = "wf_statoff", type = "boolean?", default = FALSE, position = -1)
wf_paired_stats_input <- InputParam(id = "wf_paired_stats", type = "boolean?", default = FALSE,
                                    position = -1)
wf_novelSS_input <- InputParam(id = "wf_novelSS", type = "boolean?", default = FALSE, position = -1)
wf_mil_input <- InputParam(id = "wf_mil", type = "int?", default = 50L, position = -1)
wf_mel_input <- InputParam(id = "wf_mel", type = "int?", default = 500L, position = -1)
wf_allow_clipping_input <- InputParam(id = "wf_allow_clipping", type = "boolean?", default = FALSE,
                                      position = -1)
wf_machine_mem_gb_input <- InputParam(id = "wf_machine_mem_gb", type = "int?", default = 4L,
                                      position = -1)
wf_disk_space_gb_input <- InputParam(id = "wf_disk_space_gb", type = "int?", default = 20L,
                                     position = -1)

## Pre step cwlProcess
pre_bam_input <- InputParam(id = "pre_bam", type = "File", position = -1)
pre_bam_id_input <- InputParam(id = "pre_bam_id", type = "string", position = -1)
pre_gtf_input <- InputParam(id = "pre_gtf", type = "File", position = -1)
pre_pairedORsingle_input <- InputParam(id = "pre_pairedORsingle", type = "string", position = -1)
pre_readLength_input <- InputParam(id = "pre_readLength", type = "int", position = -1)
pre_out_dir_input <- InputParam(id = "pre_out_dir", type = "string", position = -1)
pre_lib_type_input <- InputParam(id = "pre_lib_type", type = "string", position = -1)
pre_variable_readLength_input <- InputParam(id = "pre_variable_readLength", type = "boolean",
                                            position = -1)
pre_anchorLength_input <- InputParam(id = "pre_anchorLength", type = "int?", position = -1)
pre_novelSS_input <- InputParam(id = "pre_novelSS", type = "boolean", position = -1)
pre_mil_input <- InputParam(id = "pre_mil", type = "int", position = -1)
pre_mel_input <- InputParam(id = "pre_mel", type = "int", position = -1)
pre_allow_clipping_input <- InputParam(id = "pre_allow_clipping", type = "boolean", position = -1)
pre_machine_mem_gb_input <- InputParam(id = "pre_machine_mem_gb", type = "int", position = -1)
pre_disk_space_gb_input <- InputParam(id = "pre_disk_space_gb", type = "int", position = -1)
pre_script_string <- paste(sep = "\n",
"${",
"var variable_readLength_opt = inputs.pre_variable_readLength ? '--variable-read-length' : ''",
"var anchorLength_opt = inputs.pre_anchorLength != null ? '--anchorLength' : ''",
"var anchorLength_string = inputs.pre_anchorLength != null ? inputs.pre_anchorLength : ''",
"var novelSS_opt = inputs.pre_novelSS ? '--novelSS' : ''",
"var mil_opt = inputs.pre_novelSS ? '--mil' : ''",
"var mil_val = inputs.pre_novelSS ? inputs.pre_mil : ''",
"var mel_opt = inputs.pre_novelSS ? '--mel' : ''",
"var mel_val = inputs.pre_novelSS ? inputs.pre_mel : ''",
"var allow_clipping_opt = inputs.pre_allow_clipping ? '--allow-clipping' : ''",
"var script = '#!/bin/bash\\n'",
"script += 'echo ' + inputs.pre_bam.path + ' > prep.txt\\n'",
"script += 'python /rmats/rmats.py --b1 prep.txt --gtf ' + inputs.pre_gtf.path + ' -t ' + inputs.pre_pairedORsingle + ' --readLength ' + inputs.pre_readLength + ' --nthread 1 --od ' + inputs.pre_out_dir + ' --tmp tmp_output_prep_' + inputs.pre_bam_id + ' --task prep --libType ' + inputs.pre_lib_type + ' ' + variable_readLength_opt + ' ' + anchorLength_opt + ' ' + anchorLength_string + ' ' + novelSS_opt + ' ' + mil_opt + ' ' + mil_val + ' ' + mel_opt + ' ' + mel_val + ' ' + allow_clipping_opt + '\\n'",
"script += 'mkdir outfd\\n'",
"script += 'python /rmats/cp_with_prefix.py prep_' + inputs.pre_bam_id + '_ outfd tmp_output_prep_' + inputs.pre_bam_id + '/*.rmats\\n'",
"script += 'cp tmp_output_prep_' + inputs.pre_bam_id + '/*read_outcomes_by_bam.txt prep_' + inputs.pre_bam_id + '_read_outcomes_by_bam.txt\\n'",
"return(script)}")

pre_script_dirent <- Dirent(entryname = "script.sh", entry = pre_script_string)
pre_docker_req <- requireDocker(rmats_docker)
pre_js_req <- requireJS()
pre_init_work_dir_req <- requireInitialWorkDir(listing = list(pre_script_dirent))
pre_resource_req <- requireResource(coresMin = 1L, ramMin = "$(inputs.pre_machine_mem_gb * 1024)",
                                    outdirMin="$(inputs.pre_disk_space_gb * 1024)")
pre_out_rmats_output <- OutputParam(id = "pre_out_rmats", type = "File[]", glob = "outfd/*.rmats")
pre_read_outcome_output <- OutputParam(id = "pre_read_outcome", type = "File",
                                       glob = "$('prep_' + inputs.pre_bam_id + '_read_outcomes_by_bam.txt')")
pre_bam_name_output <- OutputParam(id = "pre_bam_name", type = "string",
                                   outputEval = "$(inputs.pre_bam.path.split('/').pop())")
rmats_pre <- cwlProcess(baseCommand = "bash script.sh",
                        requirements = list(pre_docker_req, pre_js_req, pre_init_work_dir_req,
                                            pre_resource_req),
                        inputs = InputParamList(pre_bam_input, pre_bam_id_input, pre_gtf_input,
                                                pre_pairedORsingle_input, pre_readLength_input,
                                                pre_out_dir_input, pre_lib_type_input,
                                                pre_variable_readLength_input,
                                                pre_anchorLength_input, pre_novelSS_input,
                                                pre_mil_input, pre_mel_input,
                                                pre_allow_clipping_input, pre_machine_mem_gb_input,
                                                pre_disk_space_gb_input),
                        outputs = OutputParamList(pre_out_rmats_output, pre_read_outcome_output,
                                                  pre_bam_name_output))

## Expression tool steps to convert Files to locations.
## This avoids loading all the Files to the disk of a single worker machine.
## .location must be used rather than .path for ExpressionTool since it should not
## have access to a real file path (although Cavatica seems to load the file to disk anyway).
exp_file_to_loc_file_input <- InputParam(id = "exp_file_to_loc_file", type = "File", position = -1)
exp_file_to_loc_js_req <- requireJS()
exp_file_to_loc_js <- "${return({'exp_file_to_loc_loc': inputs.exp_file_to_loc_file.location})}"
exp_file_to_loc_loc_output <- OutputParam(id = "exp_file_to_loc_loc", type = "string")
exp_file_to_loc <- cwlProcess(cwlClass = "ExpressionTool",
                              requirements = list(exp_file_to_loc_js_req),
                              inputs = InputParamList(exp_file_to_loc_file_input),
                              outputs = OutputParamList(exp_file_to_loc_loc_output),
                              expression = exp_file_to_loc_js)

step_exp_file_to_loc_g1 <- cwlStep(id = "step_exp_file_to_loc_g1",
                                   run = exp_file_to_loc,
                                   In = list(exp_file_to_loc_file = "wf_bam_g1"),
                                   scatter = list("exp_file_to_loc_file"),
                                   scatterMethod = "dotproduct")

step_exp_file_to_loc_g2 <- cwlStep(id = "step_exp_file_to_loc_g2",
                                   run = exp_file_to_loc,
                                   In = list(exp_file_to_loc_file = "wf_bam_g2"),
                                   scatter = list("exp_file_to_loc_file"),
                                   scatterMethod = "dotproduct")

## Expression tool steps to generate bam_ids
exp_bam_id_bams_input <- InputParam(id = "exp_bam_id_bams", type = "string[]", position = -1)
exp_bam_id_prefix_input <- InputParam(id = "exp_bam_id_prefix", type = "string", position = -1)
exp_bam_id_js_req <- requireJS()
exp_bam_id_ids_js <- paste(sep = "\n",
"${",
"var id_strings = new Array(inputs.exp_bam_id_bams.length)",
"for (var i = 0; i < id_strings.length; i++) {",
"  id_strings[i] = inputs.exp_bam_id_prefix + i",
"}",
"return({'exp_bam_id_ids': id_strings})}")
exp_bam_id_ids_output <- OutputParam(id = "exp_bam_id_ids", type = "string[]")
exp_bam_id <- cwlProcess(cwlClass = "ExpressionTool",
                         requirements = list(exp_bam_id_js_req),
                         inputs = InputParamList(exp_bam_id_bams_input, exp_bam_id_prefix_input),
                         outputs = OutputParamList(exp_bam_id_ids_output),
                         expression = exp_bam_id_ids_js)

step_exp_pre_g1 <- cwlStep(id = "step_exp_pre_g1", run = exp_bam_id,
                           In = list(exp_bam_id_bams = "step_exp_file_to_loc_g1/exp_file_to_loc_loc",
                                     exp_bam_id_prefix = list(valueFrom = "g1_")))
step_exp_pre_g2 <- cwlStep(id = "step_exp_pre_g2", run = exp_bam_id,
                           In = list(exp_bam_id_bams = "step_exp_file_to_loc_g2/exp_file_to_loc_loc",
                                     exp_bam_id_prefix = list(valueFrom = "g2_")))

## Pre step scatter over bam_g1
step_pre_g1 <- cwlStep(id = "step_pre_g1",
                       run = rmats_pre,
                       In = list(pre_bam = "wf_bam_g1",
                                 pre_bam_id = "step_exp_pre_g1/exp_bam_id_ids",
                                 pre_gtf = "wf_gtf",
                                 pre_pairedORsingle = "wf_pairedORsingle",
                                 pre_readLength = "wf_readLength",
                                 pre_out_dir = "wf_out_dir",
                                 pre_lib_type = "wf_lib_type",
                                 pre_variable_readLength = "wf_variable_readLength",
                                 pre_anchorLength = "wf_anchorLength",
                                 pre_novelSS = "wf_novelSS",
                                 pre_mil = "wf_mil",
                                 pre_mel = "wf_mel",
                                 pre_allow_clipping = "wf_allow_clipping",
                                 pre_machine_mem_gb = "wf_machine_mem_gb",
                                 pre_disk_space_gb = "wf_disk_space_gb"),
                       scatter = list("pre_bam", "pre_bam_id"),
                       scatterMethod = "dotproduct")

## Pre step scatter over bam_g2
step_pre_g2 <- cwlStep(id = "step_pre_g2",
                       run = rmats_pre,
                       In = list(pre_bam = "wf_bam_g2",
                                 pre_bam_id = "step_exp_pre_g2/exp_bam_id_ids",
                                 pre_gtf = "wf_gtf",
                                 pre_pairedORsingle = "wf_pairedORsingle",
                                 pre_readLength = "wf_readLength",
                                 pre_out_dir = "wf_out_dir",
                                 pre_lib_type = "wf_lib_type",
                                 pre_variable_readLength = "wf_variable_readLength",
                                 pre_anchorLength = "wf_anchorLength",
                                 pre_novelSS = "wf_novelSS",
                                 pre_mil = "wf_mil",
                                 pre_mel = "wf_mel",
                                 pre_allow_clipping = "wf_allow_clipping",
                                 pre_machine_mem_gb = "wf_machine_mem_gb",
                                 pre_disk_space_gb = "wf_disk_space_gb"),
                       scatter = list("pre_bam", "pre_bam_id"),
                       scatterMethod = "dotproduct")

## Post step cwlProcess
post_bam_name_g1_input <- InputParam(id = "post_bam_name_g1", type = "string[]", position = -1)
post_bam_name_g2_input <- InputParam(id = "post_bam_name_g2", type = "string[]", position = -1)
post_gtf_input <- InputParam(id = "post_gtf", type = "File", position = -1)
post_readLength_input <- InputParam(id = "post_readLength", type = "int", position = -1)
post_nthread_input <- InputParam(id = "post_nthread", type = "int", position = -1)
post_out_dir_input <- InputParam(id = "post_out_dir", type = "string", position = -1)
post_input_rmats1_input <- InputParam(id = "post_input_rmats1", type = get_2d_array_type('File'),
                                      position = -1)
post_input_rmats2_input <- InputParam(id = "post_input_rmats2", type = get_2d_array_type('File'),
                                      position = -1)
post_anchorLength_input <- InputParam(id = "post_anchorLength", type = "int?", position = -1)
post_tstat_input <- InputParam(id = "post_tstat", type = "int", position = -1)
post_cstat_input <- InputParam(id = "post_cstat", type = "float", position = -1)
post_statoff_input <- InputParam(id = "post_statoff", type = "boolean", position = -1)
post_paired_stats_input <- InputParam(id = "post_paired_stats", type = "boolean", position = -1)
post_novelSS_input <- InputParam(id = "post_novelSS", type = "boolean", position = -1)
post_mil_input <- InputParam(id = "post_mil", type = "int", position = -1)
post_mel_input <- InputParam(id = "post_mel", type = "int", position = -1)
post_machine_mem_gb_input <- InputParam(id = "post_machine_mem_gb", type = "int", position = -1)
post_disk_space_gb_input <- InputParam(id = "post_disk_space_gb", type = "int", position = -1)

post_script_string <- paste(sep = "\n",
"${",
"var anchorLength_opt = inputs.pre_anchorLength != null ? '--anchorLength' : ''",
"var anchorLength_string = inputs.pre_anchorLength != null ? inputs.pre_anchorLength : ''",
"var cstat_opt = inputs.post_paired_stats ? '' : '--cstat'",
"var cstat_val = inputs.post_paired_stats ? '' : inputs.post_cstat",
"var statoff_opt = inputs.post_statoff ? '--statoff' : ''",
"var paired_stats_opt = inputs.post_paired_stats ? '--paired-stats' : ''",
"var novelSS_opt = inputs.pre_novelSS ? '--novelSS' : ''",
"var mil_opt = inputs.pre_novelSS ? '--mil' : ''",
"var mil_val = inputs.pre_novelSS ? inputs.pre_mil : ''",
"var mel_opt = inputs.pre_novelSS ? '--mel' : ''",
"var mel_val = inputs.pre_novelSS ? inputs.pre_mel : ''",
"var rmats1 = new Array()",
"for (var i = 0; i < inputs.post_input_rmats1.length; i++) {",
"  for (var j = 0; j < inputs.post_input_rmats1[i].length; j++) {",
"    rmats1.push(inputs.post_input_rmats1[i][j])",
"  }",
"}",
"var rmats2 = new Array()",
"for (var i = 0; i < inputs.post_input_rmats2.length; i++) {",
"  for (var j = 0; j < inputs.post_input_rmats2[i].length; j++) {",
"    rmats2.push(inputs.post_input_rmats2[i][j])",
"  }",
"}",
"var rmats = new Array()",
"for (var i = 0; i < rmats1.length; i++) {",
"  rmats.push(rmats1[i])",
"}",
"for (var i = 0; i < rmats2.length; i++) {",
"  rmats.push(rmats2[i])",
"}",
"var script = '#!/bin/bash\\n'",
"script += 'mkdir fd_rmats\\n'",
"script += 'for file in'",
"for (var i = 0; i < rmats.length; i++) {",
"  script += ' ' + rmats[i].path",
"}",
## sed argument is single quoted and includes backslashes
## sed 's/.*\///g' -> sed \\'s/.*\\\\///g\\'
"script += '; do fn=`basename $file`; sed \\'s/.*\\\\///g\\' $file > fd_rmats/$fn; done\\n'",
## echo ${sep="," bam_name_g1} > bam_g1.txt
"script += 'echo '",
"for (var i = 0; i < inputs.post_bam_name_g1.length; i++) {",
"  script += inputs.post_bam_name_g1[i]",
"  if (i != (inputs.post_bam_name_g1.length - 1)) {",
"    script += ','",
"  }",
"}",
"script += ' > bam_g1.txt\\n'",
## echo ${sep="," bam_name_g2} > bam_g2.txt
"script += 'echo '",
"for (var i = 0; i < inputs.post_bam_name_g2.length; i++) {",
"  script += inputs.post_bam_name_g2[i]",
"  if (i != (inputs.post_bam_name_g2.length - 1)) {",
"    script += ','",
"  }",
"}",
"script += ' > bam_g2.txt\\n'",
"script += 'python /rmats/rmats.py --b1 bam_g1.txt --b2 bam_g2.txt --gtf ' + inputs.post_gtf.path + ' --readLength ' + inputs.post_readLength + ' --nthread ' + inputs.post_nthread + ' --od ' + inputs.post_out_dir + ' --tmp fd_rmats --task post ' + anchorLength_opt + ' ' + anchorLength_string + ' --tstat ' + inputs.post_tstat + ' ' + cstat_opt + ' ' + cstat_val + ' ' + statoff_opt + ' ' + paired_stats_opt + ' ' + novelSS_opt + ' ' + mil_opt + ' ' + mil_val + ' ' + mel_opt + ' ' + mel_val + '\\n'",
"script += 'tar czf ' + inputs.post_out_dir + '.tar.gz ' + inputs.post_out_dir + '\\n'",
"return(script)}")

post_script_dirent <- Dirent(entryname = "script.sh", entry = post_script_string)
post_docker_req <- requireDocker(rmats_docker)
post_js_req <- requireJS()
post_init_work_dir_req <- requireInitialWorkDir(listing = list(post_script_dirent))
post_resource_req <- requireResource(coresMin = "$(inputs.post_nthread)",
                                     ramMin = "$(inputs.post_machine_mem_gb * 1024)",
                                     outdirMin = "$(inputs.post_disk_space_gb * 1024)")
post_out_tar_output <- OutputParam(id = "post_out_tar", type = "File",
                                   glob = "$(inputs.post_out_dir + '.tar.gz')")
rmats_post <- cwlProcess(baseCommand = "bash script.sh",
                         requirements = list(post_docker_req, post_js_req, post_init_work_dir_req,
                                             post_resource_req),
                         inputs = InputParamList(post_bam_name_g1_input, post_bam_name_g2_input,
                                                 post_gtf_input, post_readLength_input,
                                                 post_nthread_input, post_out_dir_input,
                                                 post_input_rmats1_input, post_input_rmats2_input,
                                                 post_anchorLength_input, post_tstat_input,
                                                 post_cstat_input, post_statoff_input,
                                                 post_paired_stats_input, post_novelSS_input,
                                                 post_mil_input, post_mel_input,
                                                 post_machine_mem_gb_input,
                                                 post_disk_space_gb_input),
                         outputs = OutputParamList(post_out_tar_output))

## Post step
step_post <- cwlStep(id = "step_post",
                     run = rmats_post,
                     In = list(post_bam_name_g1 = "step_pre_g1/pre_bam_name",
                               post_bam_name_g2 = "step_pre_g2/pre_bam_name",
                               post_gtf = "wf_gtf",
                               post_readLength = "wf_readLength",
                               post_nthread = "wf_nthread",
                               post_out_dir = "wf_out_dir",
                               post_input_rmats1 = "step_pre_g1/pre_out_rmats",
                               post_input_rmats2 = "step_pre_g2/pre_out_rmats",
                               post_anchorLength = "wf_anchorLength",
                               post_tstat = "wf_tstat",
                               post_cstat = "wf_cstat",
                               post_statoff = "wf_statoff",
                               post_paired_stats = "wf_paired_stats",
                               post_novelSS = "wf_novelSS",
                               post_mil = "wf_mil",
                               post_mel = "wf_mel",
                               post_machine_mem_gb = "wf_machine_mem_gb",
                               post_disk_space_gb = "wf_disk_space_gb"))


## Expression tool step to flatten prep step read_outcome_by_bam files
exp_read_outcome_pre_g1_input <- InputParam(id = "exp_read_outcome_pre_g1", type = "File[]",
                                            position = -1)
exp_read_outcome_pre_g2_input <- InputParam(id = "exp_read_outcome_pre_g2", type = "File[]",
                                            position = -1)
exp_read_outcome_js_req <- requireJS()
exp_read_outcome_outcomes_js <- paste(sep = "\n",
"${",
"var outcomes = new Array()",
"for (var i = 0; i < inputs.exp_read_outcome_pre_g1.length; i++) {",
"  outcomes.push(inputs.exp_read_outcome_pre_g1[i])",
"}",
"for (var i = 0; i < inputs.exp_read_outcome_pre_g2.length; i++) {",
"  outcomes.push(inputs.exp_read_outcome_pre_g2[i])",
"}",
"return({'exp_read_outcome_outcomes': outcomes})}")
exp_read_outcome_outcomes_output <- OutputParam(id = "exp_read_outcome_outcomes", type = "File[]")
exp_read_outcome <- cwlProcess(cwlClass = "ExpressionTool",
                               requirements = list(exp_read_outcome_js_req),
                               inputs = InputParamList(exp_read_outcome_pre_g1_input,
                                                       exp_read_outcome_pre_g2_input),
                               outputs = OutputParamList(exp_read_outcome_outcomes_output),
                               expression = exp_read_outcome_outcomes_js)

step_exp_read_outcome <- cwlStep(id = "step_exp_read_outcome", run = exp_read_outcome,
                                 In = list(exp_read_outcome_pre_g1 = "step_pre_g1/pre_read_outcome",
                                           exp_read_outcome_pre_g2 = "step_pre_g2/pre_read_outcome"))

## Final Workflow
wf_read_outcomes_output <- OutputParam(id = "wf_read_outcomes", type = "File[]",
                                       outputSource = "step_exp_read_outcome/exp_read_outcome_outcomes")
wf_out_tar_output <- OutputParam(id = "wf_out_tar", type = "File",
                                 outputSource = "step_post/post_out_tar")
wf_scatter_req <- requireScatter()
wf_step_input_exp_req <- requireStepInputExpression()
workflow <- cwlWorkflow(requirements = list(wf_scatter_req, wf_step_input_exp_req),
                        inputs = InputParamList(wf_bam_g1_input, wf_bam_g2_input, wf_gtf_input,
                                                wf_pairedORsingle_input, wf_readLength_input,
                                                wf_nthread_input, wf_out_dir_input,
                                                wf_lib_type_input, wf_variable_readLength_input,
                                                wf_anchorLength_input, wf_tstat_input,
                                                wf_cstat_input, wf_statoff_input,
                                                wf_paired_stats_input, wf_novelSS_input,
                                                wf_mil_input, wf_mel_input, wf_allow_clipping_input,
                                                wf_machine_mem_gb_input, wf_disk_space_gb_input),
                        outputs = OutputParamList(wf_read_outcomes_output, wf_out_tar_output))
workflow <- workflow + step_exp_file_to_loc_g1 + step_exp_file_to_loc_g2 + step_exp_pre_g1 +
    step_exp_pre_g2 + step_pre_g1 + step_pre_g2 + step_post + step_exp_read_outcome
