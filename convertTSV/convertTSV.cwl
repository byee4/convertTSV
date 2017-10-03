#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [ convertTSV.py ]

arguments: [ --input_type, "csv", --output_type, "h5",  --output, $(inputs.csvfile.nameroot).h5, --genome,  "mm10"]

inputs:

  inFile:
    type: File
    inputBinding:
      prefix: --input

  inFileType:
    type: File
    inputBinding:
      prefix: --input_type

  outFile:
    type: File
    inputBinding:
      prefix: --output

  outFileType:
    type: File
    inputBinding:
      prefix: --output_type

  genome:
    type: File
    inputBinding:
      prefix: --genome

outputs:
  outFile:
    type: File
    outputBinding:
      glob: $(inputs.csvfile.nameroot).h5

