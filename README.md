![](https://github.com/dnanexus/dxCompiler/workflows/Unit%20Tests/badge.svg)
![](https://github.com/dnanexus/dxCompiler/workflows/WDL%20Integration%20Tests/badge.svg)

dxCompiler takes a pipeline written in the
[Workflow Description Language (WDL)](http://www.openwdl.org/) or [Common Workflow Language](https://www.commonwl.org/v1.2) and compiles it to an equivalent workflow on the DNAnexus platform. WDL draft-2 and versions 1.0 and 1.1 are fully supported, while WDL 2.0 (aka 'development') and CWL 1.2 support are under active development.

## Setup

Prerequisites:
* [DNAnexus platform](https://platform.dnanexus.com) account
* [dx-toolkit](https://documentation.dnanexus.com/downloads)
* java 8+
* python 2.7/3.x.

Make sure you've installed the dx-toolkit CLI, and initialized it with `dx login`. Download the latest compiler jar file from the [releases](https://github.com/dnanexus/dxCompiler/releases) page.

## Example workflow

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
    samtools index ${bam}
    mkdir slices/
    for i in `seq ${num_chrom}`; do
        samtools view -b ${bam} -o slices/$i.bam $i
    done
    >>>
    runtime {
        docker: "quay.io/ucsc_cgl/samtools"
    }
    output {
        File bai = "${bam}.bai"
        Array[File] slices = glob("slices/*.bam")
    }
}

task count_bam {
    input {
        File bam
    }

    command <<<
        samtools view -c ${bam}
    >>>
    runtime {
        docker: "quay.io/ucsc_cgl/samtools"
    }
    output {
        Int count = read_int(stdout())
    }
}
```

From the command line, we can compile the workflow to the DNAnexus platform using the dxCompiler jar file.

```
$ java -jar dxCompiler.jar compile bam_chrom_counter.wdl -project project-xxxx
```

This compiles the source WDL file to several platform objects.
- A workflow `bam_chrom_counter`
- Two applets that can be called independently: `slice_bam`, and `count_bam`
- A few auxiliary applets that process workflow inputs, outputs, and launch the scatter.

These objects are all created in the current `dx` project and folder. The generated workflow can be run using `dx run`. For example:

```
dx run bam_chrom_counter -i0.file=file-xxxx
```

The compiled workflow can be executed via the DNAnexus command line client or web UI:

![this](doc/bam_chrom_counter.png)

## Strict syntax

dxCompiler uses [wdlTools](https://github.com/dnanexus-rnd/wdlTools), a parser that adheres strictly to the WDL specifications. Most of the problematic automatic type conversions that are allowed by some other WDL runtime engines are not allowed by dxCompiler. Please use the command line tools in wdlTools (e.g. `check` and `lint`) to validate your WDL files before trying to compile them with dxCompiler.

## Limitations

* Calls with missing arguments have limited support.
* WDL Workflows with forward references are not yet supported.
* The new `Directory` type in WDL development/2.0 is not yet supported.

## Additional information

- [Advanced options](doc/ExpertOptions.md) explains additional compiler options
- [Internals](doc/Internals.md) describes current compiler structure (_work in progress_)
- [Tips](doc/Tips.md) examples for how to write good WDL code
- A high-level [list of changes](doc/WdlVersionChanges.md) between WDL draft-2 and version 1.0

## Contributing to dxCompiler

See the [development guide](doc/Developing.md#/) for more information on how to set up your development environment to contribute to dxCompiler and how to test your updates.

## Contributions

This software is a community effort! You can browse any of the contributions, that are not a part of dxCompiler main source codebase, below in our [contrib](contrib) folder, and add your own (see [Contributing to dxCompiler](#Contributing-to-dxCompiler)).

## Issues and feature requests

[Let us know](https://github.com/dnanexus/dxCompiler/issues) if you would like to contribute, request a feature, or report a bug.
