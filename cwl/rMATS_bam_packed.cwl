class: Workflow
cwlVersion: v1.2
inputs:
  wf_allow_clipping:
    default: false
    type: boolean?
  wf_anchorLength:
    type: int?
  wf_bam_g1:
    type: File[]
  wf_bam_g2:
    default: []
    type: File[]?
  wf_cstat:
    default: '0.0001'
    type: string?
  wf_darts_cutoff:
    default: '0.05'
    type: string?
  wf_darts_model:
    default: false
    type: boolean?
  wf_disk_space_gb:
    default: 20
    type: int?
  wf_gtf:
    type: File
  wf_individual_counts:
    default: false
    type: boolean?
  wf_is_single_end:
    default: false
    type: boolean?
  wf_lib_type:
    default: fr-unstranded
    type: string?
  wf_machine_mem_gb:
    default: 4
    type: int?
  wf_mel:
    default: 500
    type: int?
  wf_mil:
    default: 50
    type: int?
  wf_novelSS:
    default: false
    type: boolean?
  wf_nthread:
    default: 1
    type: int?
  wf_out_dir:
    type: string
  wf_paired_stats:
    default: false
    type: boolean?
  wf_readLength:
    type: int
  wf_statoff:
    default: false
    type: boolean?
  wf_tstat:
    default: 1
    type: int?
  wf_variable_read_length:
    default: false
    type: boolean?
outputs:
  wf_out_tar:
    outputSource: step_post/post_out_tar
    type: File
  wf_read_outcomes:
    outputSource: step_exp_read_outcome/exp_read_outcome_outcomes
    type: File[]
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
steps:
- id: step_exp_file_to_loc_g1
  in:
    exp_file_to_loc_file: wf_bam_g1
  out:
  - exp_file_to_loc_loc
  run:
    class: ExpressionTool
    cwlVersion: v1.2
    expression: '${return({''exp_file_to_loc_loc'': inputs.exp_file_to_loc_file.location})}'
    inputs:
      exp_file_to_loc_file:
        type: File
    outputs:
      exp_file_to_loc_loc:
        type: string
    requirements:
    - class: InlineJavascriptRequirement
  scatter:
  - exp_file_to_loc_file
  scatterMethod: dotproduct
- id: step_exp_file_to_loc_g2
  in:
    exp_file_to_loc_file: wf_bam_g2
  out:
  - exp_file_to_loc_loc
  run:
    class: ExpressionTool
    cwlVersion: v1.2
    expression: '${return({''exp_file_to_loc_loc'': inputs.exp_file_to_loc_file.location})}'
    inputs:
      exp_file_to_loc_file:
        type: File
    outputs:
      exp_file_to_loc_loc:
        type: string
    requirements:
    - class: InlineJavascriptRequirement
  scatter:
  - exp_file_to_loc_file
  scatterMethod: dotproduct
- id: step_exp_prep_g1
  in:
    exp_bam_id_bams: step_exp_file_to_loc_g1/exp_file_to_loc_loc
    exp_bam_id_prefix:
      id: exp_bam_id_prefix
      valueFrom: g1_
  out:
  - exp_bam_id_ids
  run:
    class: ExpressionTool
    cwlVersion: v1.2
    expression: "${\nvar id_strings = new Array(inputs.exp_bam_id_bams.length)\nfor\
      \ (var i = 0; i < id_strings.length; i++) {\n  id_strings[i] = inputs.exp_bam_id_prefix\
      \ + i\n}\nreturn({'exp_bam_id_ids': id_strings})}"
    inputs:
      exp_bam_id_bams:
        type: string[]
      exp_bam_id_prefix:
        type: string
    outputs:
      exp_bam_id_ids:
        type: string[]
    requirements:
    - class: InlineJavascriptRequirement
- id: step_exp_prep_g2
  in:
    exp_bam_id_bams: step_exp_file_to_loc_g2/exp_file_to_loc_loc
    exp_bam_id_prefix:
      id: exp_bam_id_prefix
      valueFrom: g2_
  out:
  - exp_bam_id_ids
  run:
    class: ExpressionTool
    cwlVersion: v1.2
    expression: "${\nvar id_strings = new Array(inputs.exp_bam_id_bams.length)\nfor\
      \ (var i = 0; i < id_strings.length; i++) {\n  id_strings[i] = inputs.exp_bam_id_prefix\
      \ + i\n}\nreturn({'exp_bam_id_ids': id_strings})}"
    inputs:
      exp_bam_id_bams:
        type: string[]
      exp_bam_id_prefix:
        type: string
    outputs:
      exp_bam_id_ids:
        type: string[]
    requirements:
    - class: InlineJavascriptRequirement
- id: step_prep_g1
  in:
    prep_allow_clipping: wf_allow_clipping
    prep_anchorLength: wf_anchorLength
    prep_bam: wf_bam_g1
    prep_bam_id: step_exp_prep_g1/exp_bam_id_ids
    prep_disk_space_gb: wf_disk_space_gb
    prep_gtf: wf_gtf
    prep_is_single_end: wf_is_single_end
    prep_lib_type: wf_lib_type
    prep_machine_mem_gb: wf_machine_mem_gb
    prep_mel: wf_mel
    prep_mil: wf_mil
    prep_novelSS: wf_novelSS
    prep_out_dir: wf_out_dir
    prep_readLength: wf_readLength
    prep_variable_read_length: wf_variable_read_length
  out:
  - prep_out_rmats
  - prep_read_outcome
  - prep_bam_name
  run:
    baseCommand: bash script.sh
    class: CommandLineTool
    cwlVersion: v1.2
    inputs:
      prep_allow_clipping:
        type: boolean
      prep_anchorLength:
        type: int?
      prep_bam:
        type: File
      prep_bam_id:
        type: string
      prep_disk_space_gb:
        type: int
      prep_gtf:
        type: File
      prep_is_single_end:
        type: boolean
      prep_lib_type:
        type: string
      prep_machine_mem_gb:
        type: int
      prep_mel:
        type: int
      prep_mil:
        type: int
      prep_novelSS:
        type: boolean
      prep_out_dir:
        type: string
      prep_readLength:
        type: int
      prep_variable_read_length:
        type: boolean
    outputs:
      prep_bam_name:
        outputBinding:
          outputEval: $(inputs.prep_bam.path.split('/').pop())
        type: string
      prep_out_rmats:
        outputBinding:
          glob: outfd/*.rmats
        type: File[]
      prep_read_outcome:
        outputBinding:
          glob: $('prep_' + inputs.prep_bam_id + '_read_outcomes_by_bam.txt')
        type: File
    requirements:
    - class: DockerRequirement
      dockerPull: xinglab/rmats:v4.3.0
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: '${

          var read_type_value = inputs.prep_is_single_end ? ''single'' : ''paired''

          var variable_read_length_opt = inputs.prep_variable_read_length ? ''--variable-read-length''
          : ''''

          var anchorLength_opt = inputs.prep_anchorLength != null ? ''--anchorLength''
          : ''''

          var anchorLength_string = inputs.prep_anchorLength != null ? inputs.prep_anchorLength
          : ''''

          var novelSS_opt = inputs.prep_novelSS ? ''--novelSS'' : ''''

          var mil_opt = inputs.prep_novelSS ? ''--mil'' : ''''

          var mil_val = inputs.prep_novelSS ? inputs.prep_mil : ''''

          var mel_opt = inputs.prep_novelSS ? ''--mel'' : ''''

          var mel_val = inputs.prep_novelSS ? inputs.prep_mel : ''''

          var allow_clipping_opt = inputs.prep_allow_clipping ? ''--allow-clipping''
          : ''''

          var script = ''#!/bin/bash\n''

          script += ''echo '' + inputs.prep_bam.path + '' > prep.txt\n''

          script += ''python /rmats/rmats.py --b1 prep.txt --gtf '' + inputs.prep_gtf.path
          + '' -t '' + read_type_value + '' --readLength '' + inputs.prep_readLength
          + '' --nthread 1 --od '' + inputs.prep_out_dir + '' --tmp tmp_output_prep_''
          + inputs.prep_bam_id + '' --task prep --libType '' + inputs.prep_lib_type
          + '' '' + variable_read_length_opt + '' '' + anchorLength_opt + '' '' +
          anchorLength_string + '' '' + novelSS_opt + '' '' + mil_opt + '' '' + mil_val
          + '' '' + mel_opt + '' '' + mel_val + '' '' + allow_clipping_opt + ''\n''

          script += ''mkdir outfd\n''

          script += ''python /rmats/cp_with_prefix.py prep_'' + inputs.prep_bam_id
          + ''_ outfd tmp_output_prep_'' + inputs.prep_bam_id + ''/*.rmats\n''

          script += ''cp tmp_output_prep_'' + inputs.prep_bam_id + ''/*read_outcomes_by_bam.txt
          prep_'' + inputs.prep_bam_id + ''_read_outcomes_by_bam.txt\n''

          return(script)}'
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: 1
      outdirMin: $(inputs.prep_disk_space_gb * 1024)
      ramMin: $(inputs.prep_machine_mem_gb * 1024)
  scatter:
  - prep_bam
  - prep_bam_id
  scatterMethod: dotproduct
- id: step_prep_g2
  in:
    prep_allow_clipping: wf_allow_clipping
    prep_anchorLength: wf_anchorLength
    prep_bam: wf_bam_g2
    prep_bam_id: step_exp_prep_g2/exp_bam_id_ids
    prep_disk_space_gb: wf_disk_space_gb
    prep_gtf: wf_gtf
    prep_is_single_end: wf_is_single_end
    prep_lib_type: wf_lib_type
    prep_machine_mem_gb: wf_machine_mem_gb
    prep_mel: wf_mel
    prep_mil: wf_mil
    prep_novelSS: wf_novelSS
    prep_out_dir: wf_out_dir
    prep_readLength: wf_readLength
    prep_variable_read_length: wf_variable_read_length
  out:
  - prep_out_rmats
  - prep_read_outcome
  - prep_bam_name
  run:
    baseCommand: bash script.sh
    class: CommandLineTool
    cwlVersion: v1.2
    inputs:
      prep_allow_clipping:
        type: boolean
      prep_anchorLength:
        type: int?
      prep_bam:
        type: File
      prep_bam_id:
        type: string
      prep_disk_space_gb:
        type: int
      prep_gtf:
        type: File
      prep_is_single_end:
        type: boolean
      prep_lib_type:
        type: string
      prep_machine_mem_gb:
        type: int
      prep_mel:
        type: int
      prep_mil:
        type: int
      prep_novelSS:
        type: boolean
      prep_out_dir:
        type: string
      prep_readLength:
        type: int
      prep_variable_read_length:
        type: boolean
    outputs:
      prep_bam_name:
        outputBinding:
          outputEval: $(inputs.prep_bam.path.split('/').pop())
        type: string
      prep_out_rmats:
        outputBinding:
          glob: outfd/*.rmats
        type: File[]
      prep_read_outcome:
        outputBinding:
          glob: $('prep_' + inputs.prep_bam_id + '_read_outcomes_by_bam.txt')
        type: File
    requirements:
    - class: DockerRequirement
      dockerPull: xinglab/rmats:v4.3.0
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: '${

          var read_type_value = inputs.prep_is_single_end ? ''single'' : ''paired''

          var variable_read_length_opt = inputs.prep_variable_read_length ? ''--variable-read-length''
          : ''''

          var anchorLength_opt = inputs.prep_anchorLength != null ? ''--anchorLength''
          : ''''

          var anchorLength_string = inputs.prep_anchorLength != null ? inputs.prep_anchorLength
          : ''''

          var novelSS_opt = inputs.prep_novelSS ? ''--novelSS'' : ''''

          var mil_opt = inputs.prep_novelSS ? ''--mil'' : ''''

          var mil_val = inputs.prep_novelSS ? inputs.prep_mil : ''''

          var mel_opt = inputs.prep_novelSS ? ''--mel'' : ''''

          var mel_val = inputs.prep_novelSS ? inputs.prep_mel : ''''

          var allow_clipping_opt = inputs.prep_allow_clipping ? ''--allow-clipping''
          : ''''

          var script = ''#!/bin/bash\n''

          script += ''echo '' + inputs.prep_bam.path + '' > prep.txt\n''

          script += ''python /rmats/rmats.py --b1 prep.txt --gtf '' + inputs.prep_gtf.path
          + '' -t '' + read_type_value + '' --readLength '' + inputs.prep_readLength
          + '' --nthread 1 --od '' + inputs.prep_out_dir + '' --tmp tmp_output_prep_''
          + inputs.prep_bam_id + '' --task prep --libType '' + inputs.prep_lib_type
          + '' '' + variable_read_length_opt + '' '' + anchorLength_opt + '' '' +
          anchorLength_string + '' '' + novelSS_opt + '' '' + mil_opt + '' '' + mil_val
          + '' '' + mel_opt + '' '' + mel_val + '' '' + allow_clipping_opt + ''\n''

          script += ''mkdir outfd\n''

          script += ''python /rmats/cp_with_prefix.py prep_'' + inputs.prep_bam_id
          + ''_ outfd tmp_output_prep_'' + inputs.prep_bam_id + ''/*.rmats\n''

          script += ''cp tmp_output_prep_'' + inputs.prep_bam_id + ''/*read_outcomes_by_bam.txt
          prep_'' + inputs.prep_bam_id + ''_read_outcomes_by_bam.txt\n''

          return(script)}'
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: 1
      outdirMin: $(inputs.prep_disk_space_gb * 1024)
      ramMin: $(inputs.prep_machine_mem_gb * 1024)
  scatter:
  - prep_bam
  - prep_bam_id
  scatterMethod: dotproduct
- id: step_post
  in:
    post_anchorLength: wf_anchorLength
    post_bam_name_g1: step_prep_g1/prep_bam_name
    post_bam_name_g2: step_prep_g2/prep_bam_name
    post_cstat: wf_cstat
    post_darts_cutoff: wf_darts_cutoff
    post_darts_model: wf_darts_model
    post_disk_space_gb: wf_disk_space_gb
    post_gtf: wf_gtf
    post_individual_counts: wf_individual_counts
    post_input_rmats1: step_prep_g1/prep_out_rmats
    post_input_rmats2: step_prep_g2/prep_out_rmats
    post_machine_mem_gb: wf_machine_mem_gb
    post_mel: wf_mel
    post_mil: wf_mil
    post_novelSS: wf_novelSS
    post_nthread: wf_nthread
    post_out_dir: wf_out_dir
    post_paired_stats: wf_paired_stats
    post_readLength: wf_readLength
    post_statoff: wf_statoff
    post_tstat: wf_tstat
  out:
  - post_out_tar
  run:
    baseCommand: bash script.sh
    class: CommandLineTool
    cwlVersion: v1.2
    inputs:
      post_anchorLength:
        type: int?
      post_bam_name_g1:
        type: string[]
      post_bam_name_g2:
        type: string[]
      post_cstat:
        type: string
      post_darts_cutoff:
        type: string
      post_darts_model:
        type: boolean
      post_disk_space_gb:
        type: int
      post_gtf:
        type: File
      post_individual_counts:
        type: boolean
      post_input_rmats1:
        type:
          items:
            items: File
            type: array
          type: array
      post_input_rmats2:
        type:
          items:
            items: File
            type: array
          type: array
      post_machine_mem_gb:
        type: int
      post_mel:
        type: int
      post_mil:
        type: int
      post_novelSS:
        type: boolean
      post_nthread:
        type: int
      post_out_dir:
        type: string
      post_paired_stats:
        type: boolean
      post_readLength:
        type: int
      post_statoff:
        type: boolean
      post_tstat:
        type: int
    outputs:
      post_out_tar:
        outputBinding:
          glob: $(inputs.post_out_dir + '.tar.gz')
        type: File
    requirements:
    - class: DockerRequirement
      dockerPull: xinglab/rmats:v4.3.0
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: "${\nvar has_g2 = inputs.post_bam_name_g2.length > 0\nvar b2_opt =\
          \ has_g2 ? '--b2' : ''\nvar b2_val = has_g2 ? 'bam_g2.txt' : ''\nvar anchorLength_opt\
          \ = inputs.post_anchorLength != null ? '--anchorLength' : ''\nvar anchorLength_string\
          \ = inputs.post_anchorLength != null ? inputs.post_anchorLength : ''\nvar\
          \ is_default_stats = (!inputs.post_paired_stats) && (!inputs.post_darts_model)\n\
          var cstat_opt = is_default_stats ? '--cstat' : ''\nvar cstat_val = is_default_stats\
          \ ? inputs.post_cstat : ''\nvar statoff_opt = inputs.post_statoff ? '--statoff'\
          \ : ''\nvar paired_stats_opt = inputs.post_paired_stats ? '--paired-stats'\
          \ : ''\nvar darts_model_opt = inputs.post_darts_model ? '--darts-model'\
          \ : ''\nvar darts_cutoff_opt = inputs.post_darts_model ? '--darts-cutoff'\
          \ : ''\nvar darts_cutoff_val = inputs.post_darts_model ? inputs.post_darts_cutoff\
          \ : ''\nvar novelSS_opt = inputs.post_novelSS ? '--novelSS' : ''\nvar mil_opt\
          \ = inputs.post_novelSS ? '--mil' : ''\nvar mil_val = inputs.post_novelSS\
          \ ? inputs.post_mil : ''\nvar mel_opt = inputs.post_novelSS ? '--mel' :\
          \ ''\nvar mel_val = inputs.post_novelSS ? inputs.post_mel : ''\nvar individual_counts_opt\
          \ = inputs.post_individual_counts ? '--individual-counts' : ''\nvar rmats1\
          \ = new Array()\nfor (var i = 0; i < inputs.post_input_rmats1.length; i++)\
          \ {\n  for (var j = 0; j < inputs.post_input_rmats1[i].length; j++) {\n\
          \    rmats1.push(inputs.post_input_rmats1[i][j])\n  }\n}\nvar rmats2 = new\
          \ Array()\nfor (var i = 0; i < inputs.post_input_rmats2.length; i++) {\n\
          \  for (var j = 0; j < inputs.post_input_rmats2[i].length; j++) {\n    rmats2.push(inputs.post_input_rmats2[i][j])\n\
          \  }\n}\nvar rmats = new Array()\nfor (var i = 0; i < rmats1.length; i++)\
          \ {\n  rmats.push(rmats1[i])\n}\nfor (var i = 0; i < rmats2.length; i++)\
          \ {\n  rmats.push(rmats2[i])\n}\nvar script = '#!/bin/bash\\n'\nscript +=\
          \ 'mkdir fd_rmats\\n'\nscript += 'for file in'\nfor (var i = 0; i < rmats.length;\
          \ i++) {\n  script += ' ' + rmats[i].path\n}\nscript += '; do fn=`basename\
          \ $file`; sed \\'s/.*\\\\///g\\' $file > fd_rmats/$fn; done\\n'\nscript\
          \ += 'echo '\nfor (var i = 0; i < inputs.post_bam_name_g1.length; i++) {\n\
          \  script += inputs.post_bam_name_g1[i]\n  if (i != (inputs.post_bam_name_g1.length\
          \ - 1)) {\n    script += ','\n  }\n}\nscript += ' > bam_g1.txt\\n'\nscript\
          \ += 'echo '\nfor (var i = 0; i < inputs.post_bam_name_g2.length; i++) {\n\
          \  script += inputs.post_bam_name_g2[i]\n  if (i != (inputs.post_bam_name_g2.length\
          \ - 1)) {\n    script += ','\n  }\n}\nscript += ' > bam_g2.txt\\n'\nscript\
          \ += 'python /rmats/rmats.py --b1 bam_g1.txt ' + b2_opt + ' ' + b2_val +\
          \ ' --gtf ' + inputs.post_gtf.path + ' --readLength ' + inputs.post_readLength\
          \ + ' --nthread ' + inputs.post_nthread + ' --od ' + inputs.post_out_dir\
          \ + ' --tmp fd_rmats --task post ' + anchorLength_opt + ' ' + anchorLength_string\
          \ + ' --tstat ' + inputs.post_tstat + ' ' + cstat_opt + ' ' + cstat_val\
          \ + ' ' + statoff_opt + ' ' + paired_stats_opt + ' ' + darts_model_opt +\
          \ ' ' + darts_cutoff_opt + ' ' + darts_cutoff_val + ' ' + novelSS_opt +\
          \ ' ' + mil_opt + ' ' + mil_val + ' ' + mel_opt + ' ' + mel_val + ' ' +\
          \ individual_counts_opt + '\\n'\nscript += 'tar czf ' + inputs.post_out_dir\
          \ + '.tar.gz ' + inputs.post_out_dir + '\\n'\nreturn(script)}"
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: $(inputs.post_nthread)
      outdirMin: $(inputs.post_disk_space_gb * 1024)
      ramMin: $(inputs.post_machine_mem_gb * 1024)
- id: step_exp_read_outcome
  in:
    exp_read_outcome_prep_g1: step_prep_g1/prep_read_outcome
    exp_read_outcome_prep_g2: step_prep_g2/prep_read_outcome
  out:
  - exp_read_outcome_outcomes
  run:
    class: ExpressionTool
    cwlVersion: v1.2
    expression: "${\nvar outcomes = new Array()\nfor (var i = 0; i < inputs.exp_read_outcome_prep_g1.length;\
      \ i++) {\n  outcomes.push(inputs.exp_read_outcome_prep_g1[i])\n}\nfor (var i\
      \ = 0; i < inputs.exp_read_outcome_prep_g2.length; i++) {\n  outcomes.push(inputs.exp_read_outcome_prep_g2[i])\n\
      }\nreturn({'exp_read_outcome_outcomes': outcomes})}"
    inputs:
      exp_read_outcome_prep_g1:
        type: File[]
      exp_read_outcome_prep_g2:
        type: File[]
    outputs:
      exp_read_outcome_outcomes:
        type: File[]
    requirements:
    - class: InlineJavascriptRequirement
