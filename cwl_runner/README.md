# dx-cwl-runner

This is the DNAnexus implementation of the [cwl-runner](https://www.commonwl.org/v1.1/CommandLineTool.html#Executing_CWL_documents_as_scripts) interface.

## Pre-requisites

* `pip install poetry`
* You need `dxCompiler-X.Y.Z.jar` in the current working directory, or you need to set the `DX_COMPILER_JAR` environment variable with the path to thd dxCompiler JAR file. 

## Execution

To run a test:

```bash
$ poetry run dx-cwl-runner [--output-dir <output dir>] [--quiet] <process_file> <job_file>
```

Where `process_file` is the CWL file to execute and `job_file` is the JSON file with job inputs.

`dx-cwl-runner` has two additional arguments beyond those specifiled in the cwl-runner interface:

* `--basedir`: specifies a directory other than the current working directory where input files are located
* `--dryrun`: causes commands to be echoed to `stderr` but no actions performed
