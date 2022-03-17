cwlVersion: v1.0
id: slice_bam
class: CommandLineTool
inputs:
- id: bam
  type: File
- id: num_chrom
  default: 22
  type: int
outputs:
- id: bai
  type: File
  outputBinding:
    glob: $(inputs.bam.basename).bai
- id: slices
  type: File[]
  outputBinding:
    glob: "slices/*.bam"
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
- class: InitialWorkDirRequirement
  listing:
  - entryname: slice_bam.sh
    entry: |-
      set -ex
      samtools index $1
      mkdir slices/
      for i in `seq $2`; do
          samtools view -b $1 -o slices/$i.bam $i
      done
  - entry: $(inputs.bam)
baseCommand: ["sh", "slice_bam.sh"]
arguments:
  - position: 0
    valueFrom: $(inputs.bam.basename)
  - position: 1
    valueFrom: $(inputs.num_chrom)
hints:
- class: NetworkAccess
  networkAccess: true
- class: LoadListingRequirement
  loadListing: deep_listing