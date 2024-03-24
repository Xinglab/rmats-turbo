class: Workflow
cwlVersion: v1.0
inputs:
  wf_allow_clipping:
    default: false
    type: boolean?
  wf_anchorLength:
    type: int?
  wf_bam_g1:
    type: File[]
  wf_bam_g2:
    type: File[]
  wf_cstat:
    default: 0.0001
    type: float?
  wf_disk_space_gb:
    default: 20
    type: int?
  wf_gtf:
    type: File
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
  wf_pairedORsingle:
    default: paired
    type: string?
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
  wf_variable_readLength:
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
    cwlVersion: v1.0
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
    cwlVersion: v1.0
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
- id: step_exp_pre_g1
  in:
    exp_bam_id_bams: step_exp_file_to_loc_g1/exp_file_to_loc_loc
    exp_bam_id_prefix:
      id: exp_bam_id_prefix
      valueFrom: g1_
  out:
  - exp_bam_id_ids
  run:
    class: ExpressionTool
    cwlVersion: v1.0
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
- id: step_exp_pre_g2
  in:
    exp_bam_id_bams: step_exp_file_to_loc_g2/exp_file_to_loc_loc
    exp_bam_id_prefix:
      id: exp_bam_id_prefix
      valueFrom: g2_
  out:
  - exp_bam_id_ids
  run:
    class: ExpressionTool
    cwlVersion: v1.0
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
- id: step_pre_g1
  in:
    pre_allow_clipping: wf_allow_clipping
    pre_anchorLength: wf_anchorLength
    pre_bam: wf_bam_g1
    pre_bam_id: step_exp_pre_g1/exp_bam_id_ids
    pre_disk_space_gb: wf_disk_space_gb
    pre_gtf: wf_gtf
    pre_lib_type: wf_lib_type
    pre_machine_mem_gb: wf_machine_mem_gb
    pre_mel: wf_mel
    pre_mil: wf_mil
    pre_novelSS: wf_novelSS
    pre_out_dir: wf_out_dir
    pre_pairedORsingle: wf_pairedORsingle
    pre_readLength: wf_readLength
    pre_variable_readLength: wf_variable_readLength
  out:
  - pre_out_rmats
  - pre_read_outcome
  - pre_bam_name
  run:
    baseCommand: bash script.sh
    class: CommandLineTool
    cwlVersion: v1.0
    inputs:
      pre_allow_clipping:
        type: boolean
      pre_anchorLength:
        type: int?
      pre_bam:
        type: File
      pre_bam_id:
        type: string
      pre_disk_space_gb:
        type: int
      pre_gtf:
        type: File
      pre_lib_type:
        type: string
      pre_machine_mem_gb:
        type: int
      pre_mel:
        type: int
      pre_mil:
        type: int
      pre_novelSS:
        type: boolean
      pre_out_dir:
        type: string
      pre_pairedORsingle:
        type: string
      pre_readLength:
        type: int
      pre_variable_readLength:
        type: boolean
    outputs:
      pre_bam_name:
        outputBinding:
          outputEval: $(inputs.pre_bam.path.split('/').pop())
        type: string
      pre_out_rmats:
        outputBinding:
          glob: outfd/*.rmats
        type: File[]
      pre_read_outcome:
        outputBinding:
          glob: $('prep_' + inputs.pre_bam_id + '_read_outcomes_by_bam.txt')
        type: File
    requirements:
    - class: DockerRequirement
      dockerPull: xinglab/rmats:v4.1.2-no-entrypoint
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: '${

          var variable_readLength_opt = inputs.pre_variable_readLength ? ''--variable-read-length''
          : ''''

          var anchorLength_opt = inputs.pre_anchorLength != null ? ''--anchorLength''
          : ''''

          var anchorLength_string = inputs.pre_anchorLength != null ? inputs.pre_anchorLength
          : ''''

          var novelSS_opt = inputs.pre_novelSS ? ''--novelSS'' : ''''

          var mil_opt = inputs.pre_novelSS ? ''--mil'' : ''''

          var mil_val = inputs.pre_novelSS ? inputs.pre_mil : ''''

          var mel_opt = inputs.pre_novelSS ? ''--mel'' : ''''

          var mel_val = inputs.pre_novelSS ? inputs.pre_mel : ''''

          var allow_clipping_opt = inputs.pre_allow_clipping ? ''--allow-clipping''
          : ''''

          var script = ''#!/bin/bash\n''

          script += ''echo '' + inputs.pre_bam.path + '' > prep.txt\n''

          script += ''python /rmats/rmats.py --b1 prep.txt --gtf '' + inputs.pre_gtf.path
          + '' -t '' + inputs.pre_pairedORsingle + '' --readLength '' + inputs.pre_readLength
          + '' --nthread 1 --od '' + inputs.pre_out_dir + '' --tmp tmp_output_prep_''
          + inputs.pre_bam_id + '' --task prep --libType '' + inputs.pre_lib_type
          + '' '' + variable_readLength_opt + '' '' + anchorLength_opt + '' '' + anchorLength_string
          + '' '' + novelSS_opt + '' '' + mil_opt + '' '' + mil_val + '' '' + mel_opt
          + '' '' + mel_val + '' '' + allow_clipping_opt + ''\n''

          script += ''mkdir outfd\n''

          script += ''python /rmats/cp_with_prefix.py prep_'' + inputs.pre_bam_id
          + ''_ outfd tmp_output_prep_'' + inputs.pre_bam_id + ''/*.rmats\n''

          script += ''cp tmp_output_prep_'' + inputs.pre_bam_id + ''/*read_outcomes_by_bam.txt
          prep_'' + inputs.pre_bam_id + ''_read_outcomes_by_bam.txt\n''

          return(script)}'
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: 1
      outdirMin: $(inputs.pre_disk_space_gb * 1024)
      ramMin: $(inputs.pre_machine_mem_gb * 1024)
  scatter:
  - pre_bam
  - pre_bam_id
  scatterMethod: dotproduct
- id: step_pre_g2
  in:
    pre_allow_clipping: wf_allow_clipping
    pre_anchorLength: wf_anchorLength
    pre_bam: wf_bam_g2
    pre_bam_id: step_exp_pre_g2/exp_bam_id_ids
    pre_disk_space_gb: wf_disk_space_gb
    pre_gtf: wf_gtf
    pre_lib_type: wf_lib_type
    pre_machine_mem_gb: wf_machine_mem_gb
    pre_mel: wf_mel
    pre_mil: wf_mil
    pre_novelSS: wf_novelSS
    pre_out_dir: wf_out_dir
    pre_pairedORsingle: wf_pairedORsingle
    pre_readLength: wf_readLength
    pre_variable_readLength: wf_variable_readLength
  out:
  - pre_out_rmats
  - pre_read_outcome
  - pre_bam_name
  run:
    baseCommand: bash script.sh
    class: CommandLineTool
    cwlVersion: v1.0
    inputs:
      pre_allow_clipping:
        type: boolean
      pre_anchorLength:
        type: int?
      pre_bam:
        type: File
      pre_bam_id:
        type: string
      pre_disk_space_gb:
        type: int
      pre_gtf:
        type: File
      pre_lib_type:
        type: string
      pre_machine_mem_gb:
        type: int
      pre_mel:
        type: int
      pre_mil:
        type: int
      pre_novelSS:
        type: boolean
      pre_out_dir:
        type: string
      pre_pairedORsingle:
        type: string
      pre_readLength:
        type: int
      pre_variable_readLength:
        type: boolean
    outputs:
      pre_bam_name:
        outputBinding:
          outputEval: $(inputs.pre_bam.path.split('/').pop())
        type: string
      pre_out_rmats:
        outputBinding:
          glob: outfd/*.rmats
        type: File[]
      pre_read_outcome:
        outputBinding:
          glob: $('prep_' + inputs.pre_bam_id + '_read_outcomes_by_bam.txt')
        type: File
    requirements:
    - class: DockerRequirement
      dockerPull: xinglab/rmats:v4.1.2-no-entrypoint
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: '${

          var variable_readLength_opt = inputs.pre_variable_readLength ? ''--variable-read-length''
          : ''''

          var anchorLength_opt = inputs.pre_anchorLength != null ? ''--anchorLength''
          : ''''

          var anchorLength_string = inputs.pre_anchorLength != null ? inputs.pre_anchorLength
          : ''''

          var novelSS_opt = inputs.pre_novelSS ? ''--novelSS'' : ''''

          var mil_opt = inputs.pre_novelSS ? ''--mil'' : ''''

          var mil_val = inputs.pre_novelSS ? inputs.pre_mil : ''''

          var mel_opt = inputs.pre_novelSS ? ''--mel'' : ''''

          var mel_val = inputs.pre_novelSS ? inputs.pre_mel : ''''

          var allow_clipping_opt = inputs.pre_allow_clipping ? ''--allow-clipping''
          : ''''

          var script = ''#!/bin/bash\n''

          script += ''echo '' + inputs.pre_bam.path + '' > prep.txt\n''

          script += ''python /rmats/rmats.py --b1 prep.txt --gtf '' + inputs.pre_gtf.path
          + '' -t '' + inputs.pre_pairedORsingle + '' --readLength '' + inputs.pre_readLength
          + '' --nthread 1 --od '' + inputs.pre_out_dir + '' --tmp tmp_output_prep_''
          + inputs.pre_bam_id + '' --task prep --libType '' + inputs.pre_lib_type
          + '' '' + variable_readLength_opt + '' '' + anchorLength_opt + '' '' + anchorLength_string
          + '' '' + novelSS_opt + '' '' + mil_opt + '' '' + mil_val + '' '' + mel_opt
          + '' '' + mel_val + '' '' + allow_clipping_opt + ''\n''

          script += ''mkdir outfd\n''

          script += ''python /rmats/cp_with_prefix.py prep_'' + inputs.pre_bam_id
          + ''_ outfd tmp_output_prep_'' + inputs.pre_bam_id + ''/*.rmats\n''

          script += ''cp tmp_output_prep_'' + inputs.pre_bam_id + ''/*read_outcomes_by_bam.txt
          prep_'' + inputs.pre_bam_id + ''_read_outcomes_by_bam.txt\n''

          return(script)}'
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: 1
      outdirMin: $(inputs.pre_disk_space_gb * 1024)
      ramMin: $(inputs.pre_machine_mem_gb * 1024)
  scatter:
  - pre_bam
  - pre_bam_id
  scatterMethod: dotproduct
- id: step_post
  in:
    post_anchorLength: wf_anchorLength
    post_bam_name_g1: step_pre_g1/pre_bam_name
    post_bam_name_g2: step_pre_g2/pre_bam_name
    post_cstat: wf_cstat
    post_disk_space_gb: wf_disk_space_gb
    post_gtf: wf_gtf
    post_input_rmats1: step_pre_g1/pre_out_rmats
    post_input_rmats2: step_pre_g2/pre_out_rmats
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
    cwlVersion: v1.0
    inputs:
      post_anchorLength:
        type: int?
      post_bam_name_g1:
        type: string[]
      post_bam_name_g2:
        type: string[]
      post_cstat:
        type: float
      post_disk_space_gb:
        type: int
      post_gtf:
        type: File
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
      dockerPull: xinglab/rmats:v4.1.2-no-entrypoint
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entry: "${\nvar anchorLength_opt = inputs.pre_anchorLength != null ? '--anchorLength'\
          \ : ''\nvar anchorLength_string = inputs.pre_anchorLength != null ? inputs.pre_anchorLength\
          \ : ''\nvar cstat_opt = inputs.post_paired_stats ? '' : '--cstat'\nvar cstat_val\
          \ = inputs.post_paired_stats ? '' : inputs.post_cstat\nvar statoff_opt =\
          \ inputs.post_statoff ? '--statoff' : ''\nvar paired_stats_opt = inputs.post_paired_stats\
          \ ? '--paired-stats' : ''\nvar novelSS_opt = inputs.pre_novelSS ? '--novelSS'\
          \ : ''\nvar mil_opt = inputs.pre_novelSS ? '--mil' : ''\nvar mil_val = inputs.pre_novelSS\
          \ ? inputs.pre_mil : ''\nvar mel_opt = inputs.pre_novelSS ? '--mel' : ''\n\
          var mel_val = inputs.pre_novelSS ? inputs.pre_mel : ''\nvar rmats1 = new\
          \ Array()\nfor (var i = 0; i < inputs.post_input_rmats1.length; i++) {\n\
          \  for (var j = 0; j < inputs.post_input_rmats1[i].length; j++) {\n    rmats1.push(inputs.post_input_rmats1[i][j])\n\
          \  }\n}\nvar rmats2 = new Array()\nfor (var i = 0; i < inputs.post_input_rmats2.length;\
          \ i++) {\n  for (var j = 0; j < inputs.post_input_rmats2[i].length; j++)\
          \ {\n    rmats2.push(inputs.post_input_rmats2[i][j])\n  }\n}\nvar rmats\
          \ = new Array()\nfor (var i = 0; i < rmats1.length; i++) {\n  rmats.push(rmats1[i])\n\
          }\nfor (var i = 0; i < rmats2.length; i++) {\n  rmats.push(rmats2[i])\n\
          }\nvar script = '#!/bin/bash\\n'\nscript += 'mkdir fd_rmats\\n'\nscript\
          \ += 'for file in'\nfor (var i = 0; i < rmats.length; i++) {\n  script +=\
          \ ' ' + rmats[i].path\n}\nscript += '; do fn=`basename $file`; sed \\'s/.*\\\
          \\///g\\' $file > fd_rmats/$fn; done\\n'\nscript += 'echo '\nfor (var i\
          \ = 0; i < inputs.post_bam_name_g1.length; i++) {\n  script += inputs.post_bam_name_g1[i]\n\
          \  if (i != (inputs.post_bam_name_g1.length - 1)) {\n    script += ','\n\
          \  }\n}\nscript += ' > bam_g1.txt\\n'\nscript += 'echo '\nfor (var i = 0;\
          \ i < inputs.post_bam_name_g2.length; i++) {\n  script += inputs.post_bam_name_g2[i]\n\
          \  if (i != (inputs.post_bam_name_g2.length - 1)) {\n    script += ','\n\
          \  }\n}\nscript += ' > bam_g2.txt\\n'\nscript += 'python /rmats/rmats.py\
          \ --b1 bam_g1.txt --b2 bam_g2.txt --gtf ' + inputs.post_gtf.path + ' --readLength\
          \ ' + inputs.post_readLength + ' --nthread ' + inputs.post_nthread + ' --od\
          \ ' + inputs.post_out_dir + ' --tmp fd_rmats --task post ' + anchorLength_opt\
          \ + ' ' + anchorLength_string + ' --tstat ' + inputs.post_tstat + ' ' +\
          \ cstat_opt + ' ' + cstat_val + ' ' + statoff_opt + ' ' + paired_stats_opt\
          \ + ' ' + novelSS_opt + ' ' + mil_opt + ' ' + mil_val + ' ' + mel_opt +\
          \ ' ' + mel_val + '\\n'\nscript += 'tar czf ' + inputs.post_out_dir + '.tar.gz\
          \ ' + inputs.post_out_dir + '\\n'\nreturn(script)}"
        entryname: script.sh
        writable: false
    - class: ResourceRequirement
      coresMin: $(inputs.post_nthread)
      outdirMin: $(inputs.post_disk_space_gb * 1024)
      ramMin: $(inputs.post_machine_mem_gb * 1024)
- id: step_exp_read_outcome
  in:
    exp_read_outcome_pre_g1: step_pre_g1/pre_read_outcome
    exp_read_outcome_pre_g2: step_pre_g2/pre_read_outcome
  out:
  - exp_read_outcome_outcomes
  run:
    class: ExpressionTool
    cwlVersion: v1.0
    expression: "${\nvar outcomes = new Array()\nfor (var i = 0; i < inputs.exp_read_outcome_pre_g1.length;\
      \ i++) {\n  outcomes.push(inputs.exp_read_outcome_pre_g1[i])\n}\nfor (var i\
      \ = 0; i < inputs.exp_read_outcome_pre_g2.length; i++) {\n  outcomes.push(inputs.exp_read_outcome_pre_g2[i])\n\
      }\nreturn({'exp_read_outcome_outcomes': outcomes})}"
    inputs:
      exp_read_outcome_pre_g1:
        type: File[]
      exp_read_outcome_pre_g2:
        type: File[]
    outputs:
      exp_read_outcome_outcomes:
        type: File[]
    requirements:
    - class: InlineJavascriptRequirement
