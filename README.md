![](https://github.com/dnanexus/dxCompiler/workflows/Unit%20Tests/badge.svg)
![](https://github.com/dnanexus/dxCompiler/workflows/WDL%20Integration%20Tests/badge.svg)

dxCompiler takes a pipeline written in the
[Workflow Description Language (WDL)](http://www.openwdl.org/) or [Common Workflow Language](https://www.commonwl.org/v1.2) and compiles it to an equivalent workflow on the DNAnexus platform. 
The following standards are fully supported:
* WDL: draft-2, 1.0, and 1.1
* CWL: 1.2 

Support for WDL 2.0 (aka 'development') is under active development and not yet production-ready. CWL 1.0 and 1.1 are not supported but can be upgraded to 1.2 to be compiled (see [Preprocess CWL workflow](#preprocess-cwl-workflow))

## Setup

To compile and run your workflow, make sure you have the following ready:
* [DNAnexus platform](https://platform.dnanexus.com) account
* [dx-toolkit](https://documentation.dnanexus.com/downloads)
  * Log in using `dx login`
  * It is recommended to pre-select the project where you want to your compiled workflows to go using `dx select`
* Java 8 or 11
* The latest dxCompiler JAR file from the [releases](https://github.com/dnanexus/dxCompiler/releases) page.
* [docker](https://docs.docker.com/get-docker/), if you want to invoke dxCompiler with the [run-dxcompiler-docker](https://github.com/dnanexus/dxCompiler/blob/main/scripts/compiler_image/run-dxcompiler-docker) script using a public `dnanexus/dxcompiler` docker container.
* Python 3.x to run the dxCompiler integration tests

To compile CWL tools/workflows, you might also need:
  * [sbpack](https://github.com/rabix/sbpack) to pack the workflow made up of multiple files into a single compound JSON document before compilation
  * [cwl-utils](https://github.com/common-workflow-language/cwl-utils) if you want to convert expression tool into commandline tool in the workflow
  * [cwl-upgrader](https://github.com/common-workflow-language/cwl-upgrader) to upgrade your workflow to version 1.2 
  * [cwltool]() to validate your workflow or test it locally 




## WDL

### Compile and run workflow

The `bam_chrom_counter` workflow is written in WDL. Task `slice_bam` splits a bam file into an array of sub-files. Task
`count_bam` counts the number of alignments on a bam file. The workflow takes an input bam file, calls `slice_bam` to split it into chromosomes, and calls `count_bam` in parallel on each chromosome. The results comprise a bam index file, and an array with the number of reads per chromosome.

```wdl
version 1.0

workflow bam_chrom_counter {
    input {
        File bam
    }

    call slice_bam {
        input : bam = bam
    }
    scatter (slice in slice_bam.slices) {
        call count_bam {
            input: bam = slice
        }
    }
    output {
        File bai = slice_bam.bai
        Array[Int] count = count_bam.count
    }
}

task slice_bam {
    input {
        File bam
        Int num_chrom = 22
    }
    command <<<
    set -ex
    samtools index ~{bam}
    mkdir slices/
    for i in `seq ~{num_chrom}`; do
        samtools view -b ~{bam} -o slices/$i.bam $i
    done
    >>>
    runtime {
        docker: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }
    output {
        File bai = "~{bam}.bai"
        Array[File] slices = glob("slices/*.bam")
    }
}

task count_bam {
    input {
        File bam
    }

    command <<<
        samtools view -c ~{bam}
    >>>
    runtime {
        docker: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }
    output {
        Int count = read_int(stdout())
    }
}
```

From the command line, we can compile the workflow to the DNAnexus platform using the dxCompiler jar file.

```
$ java -jar dxCompiler.jar compile bam_chrom_counter.wdl -project project-xxxx -folder /my/workflows/
```

This compiles the source WDL file to several platform objects.
- A workflow `bam_chrom_counter`
- Two applets that can be called independently: `slice_bam`, and `count_bam`
- A few auxiliary applets that process workflow inputs, outputs, and launch the scatter.

These objects are all created in the specified DNAnexus project and folder. The generated workflow can be run using `dx run`. For example:

```
dx run bam_chrom_counter -istage-common.bam=project-xxxx:file-yyyy
```

The compiled workflow can be executed via the DNAnexus command line client or web UI:

![this](doc/bam_chrom_counter.png)

### Validate the workflow

dxCompiler uses [wdlTools](https://github.com/dnanexus/wdlTools), a parser that adheres strictly to the WDL specifications. Most of the problematic automatic type conversions that are allowed by some other WDL runtime engines are not allowed by dxCompiler. Please use the command line tools in wdlTools (e.g. `check` and `lint`) to validate your WDL files before trying to compile them with dxCompiler.

## CWL

### Preprocess CWL workflow

dxCompiler requires the source CWL file be "packed" as a `cwl.json` file, which contains a single compound workflow with all the dependent processes included.
Before compiling your CWL workflow/tool, you might first need to process it into the packed format using the following steps:

1. Install `cwl-upgrader` and upgrade your workflow to CWL 1.2 (if it is not already): 
    ```
    $ pip install cwl-upgrader
    $ cwl-upgrader my-workflow.cwl [subworkflow1.cwl subworkflow2.cwl ...]
    ```
2. Install `sbpack` package and use the `cwlpack` command on the top-level workflow file to build the "packed" one:
    ```
    $ pip install https://github.com/rabix/sbpack.git`
    $ cwlpack --add-ids --json my-workflow.cwl > my-workflow.cwl.json
    ```
3. De-localize all local paths referenced in your CWL: if your CWL specifies a local path, e.g. a schema or a default value for a `file`-type input, you need to upload those files to a DNAnexus project and then replace the local path in your CWL with a DNAnexus URI, e.g. `dx://project-XXX:file-YYY`.

### Compile and run workflow
The same `bam_chrom_counter` workflow as the example WDL above is now written in CWL v1.0:
```
cwlVersion: v1.0
$graph:
- id: bam_chrom_counter
  class: Workflow
  requirements:
  - class: ScatterFeatureRequirement
  inputs:
  - id: bam
    type: File
  outputs:
  - id: bai
    type: File
    outputSource: slice_bam/bai
  - id: count
    type: int[]
    outputSource: count_bam/count
  steps:
  - id: slice_bam
    run: "#slice_bam"
    in:
      bam: bam
    out: [bai, slices]
  - id: count_bam
    run: "#count_bam"
    scatter: bam
    in:
      bam: slice_bam/slices
    out: [count]
- id: slice_bam
  class: CommandLineTool
  inputs:
  - id: bam
    type: File
    inputBinding:
      position: 1
  - id: num_chrom
    default: 22
    type: int
    inputBinding:
      position: 2
  outputs:
  - id: bai
    type: File
    outputBinding:
      glob: $(inputs.bam).bai
  - id: slices
    type: File[]
    outputBinding:
      glob: "*.bam"
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
  baseCommand: ["sh", "slice_bam.sh"]
- id: count_bam
  class: CommandLineTool
  requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
  inputs:
  - id: bam
    type: File
    inputBinding:
      position: 1
  baseCommand: ["samtools", "view", "-c"]
  outputs:
  - id: count
    type: int
    outputBinding:
      glob: stdout
      loadContents: true
      outputEval: "$(parseInt(self[0].contents))"
  stdout: stdout
```
Once it is upgraded and packed into [`bam_chrom_counter.cwl.json`](contrib/beginner_example/bam_chrom_counter.cwl.json) as suggested above, we can compile it to the DNAnexus platform and run it.
```
$ java -jar dxCompiler.jar compile bam_chrom_counter.cwl.json -project project-xxxx -folder /my/workflows/
$ dx run bam_chrom_counter -istage-common.bam=project-xxxx:file-yyyy
```
### Validate the workflow

dxCompiler support [CWL v1.2 standard](https://www.commonwl.org/v1.2/index.html) and compiles tools/workflows written using corresponding syntax specifications. You can use `cwltool --validate` to validate the packed CWL file you want to compile.
## Limitations

* WDL and CWL
  * Calls with missing arguments have [limited support](doc/ExpertOptions.md#task-and-workflow-inputs).
  * All task and workflow names must be unique across the entire import tree
    * For example, if `A.wdl` imports `B.wdl` and `A.wdl` defines workflow `foo`, then `B.wdl` cannot have a workflow or task named `foo`
  * Subworkflows built from higher-level workflows are not intented to be used on their own.
* WDL only
  * Workflows with forward references (i.e. a variable referenced before it is declared) are not yet supported.
  * The [alternative workflow output syntax](https://github.com/openwdl/wdl/blob/main/versions/draft-2/SPEC.md#outputs) that has been deprecated since WDL draft2 is not supported
  * The `call ... after` syntax introduced in WDL 1.1 is not yet supported

## Additional information

- [Advanced options](doc/ExpertOptions.md) explains additional compiler options
- [Internals](doc/Internals.md) describes current compiler structure (_work in progress_)
- [Tips](doc/Tips.md) examples for how to write good WDL code
- A high-level [list of changes](https://github.com/openwdl/wdl/blob/main/versions/Differences.md#draft-2-to-10) between WDL draft-2 and version 1.0

## Contributing to dxCompiler

See the [development guide](doc/Developing.md#/) for more information on how to set up your development environment to contribute to dxCompiler and how to test your updates.

## Contributions

This software is a community effort! You can browse any of the contributions, that are not a part of dxCompiler main source codebase, below in our [contrib](contrib) folder, and add your own (see [Contributing to dxCompiler](#Contributing-to-dxCompiler)).

## Issues and feature requests

[Let us know](https://github.com/dnanexus/dxCompiler/issues) if you would like to contribute, request a feature, or report a bug.
