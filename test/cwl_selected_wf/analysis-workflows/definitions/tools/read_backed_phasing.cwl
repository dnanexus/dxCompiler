#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "Read-backed phasing"
requirements:
- class: ResourceRequirement
  ramMin: 9000
  tmpdirMin: 25000
- class: DockerRequirement
  dockerPull: mgibio/gatk-cwl:3.6.0
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  reference:
    type:
    - string
    - File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-R"
      position: 1
  bam:
    type: File
    inputBinding:
      prefix: "-I"
      position: 2
    secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"}
      else {return self.basename + ".crai"}}
  vcf:
    type: File
    inputBinding:
      prefix: "-V"
      position: 3
    secondaryFiles: [".tbi"]
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T",
  "ReadBackedPhasing"]
arguments: ["-L", valueFrom: $(inputs.vcf), "-o", valueFrom: $(runtime.outdir)/phased.vcf]
outputs:
  phased_vcf:
    type: File
    outputBinding:
      glob: "phased.vcf"
