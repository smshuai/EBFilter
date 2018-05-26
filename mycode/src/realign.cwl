#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: realign.sh
inputs:
  miniBAM:
    type: File
  U1BED:
    type: File
  U1CORE:
    type: File
  ID:
    type: string
  OUT:
    type: string
outputs:
  realnBAM:
    type: File
    outputBinding:
      glob: $(inputs.out)_realign.bam
  processBAM:
    type: File
    outputBinding:
      glob: $(inputs.out)_processed.bam
  mstat:
    type: File
    outputBinding:
      glob: $(inputs.out)_mstat.csv
  umap:
    type: File
    outputBinding:
      glob: $(inputs.out)_umap.txt
  dp:
    type: File
    outputBinding:
      glob: ${inputs.out}_weighted_dp.tsv