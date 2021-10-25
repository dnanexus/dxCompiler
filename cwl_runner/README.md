# dx-cwl-runner

This is the DNAnexus implementation of the [cwl-runner](https://www.commonwl.org/v1.1/CommandLineTool.html#Executing_CWL_documents_as_scripts) interface.

## Pre-requisites

* Java 8
* Python 3.6+
* [poetry](https://python-poetry.org/) (`pip install poetry`)

## Installation

In the `cwl_runner` folder run:

```bash
$ poetry install
```

## Setup

* You must have a [DNAnexus](https://platform.dnanexus.com) account
* You must be logged into your DNAnexus account (using `dx login`)
* You must select the project where you want test data to be upload and tools/workflows to be compiled and executed (using `dx select`)
* If you want to use a folder in the project other than the root folder, you must select that folder (using `dx cd`)
* You may place `dxCompiler-X.Y.Z.jar` in the current working directory, or you may set the `DX_COMPILER_JAR` environment variable with the path to the dxCompiler JAR file. If no dxCompiler JAR can be found locally, the latest version is downloaded from [GitHub](https://github.com/dnanexus/dxCompiler/releases).

## Execution

To run a test:

```bash
$ poetry run dx-cwl-runner [--output-dir <output dir>] [--quiet] <process_file> <job_file>
```

Where `process_file` is the CWL file to execute and `job_file` is the YAML or JSON file with job inputs.

`dx-cwl-runner` has two additional options beyond those specifiled in the cwl-runner interface:

* `--basedir`: specifies a directory other than the current working directory where input files are located
* `--dryrun`: causes commands to be echoed to `stderr` but no actions performed
