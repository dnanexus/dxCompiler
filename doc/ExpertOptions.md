The reader is assumed to understand the [Workflow Description Language (WDL)](http://www.openwdl.org/) and [Common Workflow Language (CWL)](https://www.commonwl.org/v1.2), and have some experience using the [DNAnexus](http://www.dnanexus.com) platform.

dxCompiler takes a pipeline written in WDL or CWL and statically compiles it to an equivalent workflow on the DNAnexus platform. This document will use WDL examples to explain additional compiler options and features. To implement them when working with CWL workflows, please refer to [CWL v1.2.0 to WDL v1.0 mapping](CWL_v1.2.0_to_WDL_v1.md) for type and syntax equivalence between WDL and CWL.

- [Getting started](#getting-started)
  * [Compiling Workflow](#compiling-workflow)
    * [Inputs](#inputs)
    * [Defaults](#defaults)
    * [Extras](#extras)
    * [Describe WDL workflow to obtain execution tree](#describe-wdl-workflow-to-obtain-execution-tree)
- [Task and workflow inputs](#task-and-workflow-inputs)
  * [Directories](#directories)
- [Task metadata and runtime](#task-metadata-and-runtime)
  * [Task meta and parameter_meta](#task-meta-and-parameter_meta)
    * [meta section](#meta-section)
    * [parameter_meta section](#parameter_meta-section)
    * [streaming](#streaming)
  * [Task runtime and hints](#task-runtime-and-hints)
    * [Instance type](#instance-type)
    * [Additional DNAnexus-specific runtime settings](#additional-dnanexus-specific-runtime-settings)
    * [Native DNAnexus executable](#native-dnanexus-executable)
  * [Example tasks with DNAnexus-specific metadata and runtime](#example-tasks-with-dnanexus-specific-metadata-and-runtime)
- [Calling existing applets](#calling-existing-applets)
  * [Calling apps](#calling-apps)
  * [Calling apps and applets using WDL development](#calling-apps-and-applets-using-wdl-development)
  * [Overriding the native app and applet instance type](#overriding-the-native-app-and-applet-instance-type)
   * [Unsupported overrides](#unsupported-overrides)
- [Workflow metadata](#workflow-metadata)
- [Setting DNAnexus-specific attributes in extras.json](#setting-dnanexus-specific-attributes-in-extrasjson)
  * [Default and per-task attributes](#default-and-per-task-attributes)
  * [Default and per-workflow attributes](#default-and-per-workflow-attributes)
  * [Job reuse](#job-reuse)
  * [Delay workspace destruction](#delay-workspace-destruction)
- [Handling intermediate workflow outputs](#handling-intermediate-workflow-outputs)
  * [Use your own applet](#use-your-own-applet)
  * [Adding config-file based reorg applet at compilation time](#adding-config-file-based-reorg-applet-at-compilation-time)
- [Top-level calls compiled as stages](#top-level-calls-compiled-as-stages)
- [Manifests](#manifests)  
- [Docker](#docker)
  * [Setting a default docker image for all tasks](#setting-a-default-docker-image-for-all-tasks)
  * [Private registries](#private-registries)
  * [Storing a docker image as a file](#storing-a-docker-image-as-a-file)
- [Proxy configurations](#proxy-configurations)
- [Debugging an applet](#debugging-an-applet)
  * [Logging](#logging)
  * [Getting applet sources](#getting-applet-sources)
  * [Getting WDL sources](#getting-wdl-sources)
- [Recompilation](#recompilation)
- [Publishing global workflows](#publishing-global-workflows)
  * [Global workflow recommendations](#global-workflow-recommendations)
  * [Global workflow limitations](#global-workflow-limitations)


# Getting started

Prerequisites: [DNAnexus platform](https://platform.dnanexus.com) account, [dx-toolkit](https://github.com/dnanexus/dx-toolkit), java 8+, python 2.7 or 3.5+.

Make sure you've installed the dx-toolkit CLI, and initialized it with `dx login`. Download the latest dxCompiler compiler jar file from the [releases](https://github.com/dnanexus/dxCompiler/releases) page.

## Compiling Workflow

To compile a workflow:
```console
$ java -jar dxCompiler-xxx.jar compile /path/to/foo.wdl -project project-xxxx -folder /my/workflows/
```
This compiles `foo.wdl` to platform workflow `foo` in specified dx's project and folder (defaults to currently selected project and '/'). The generated workflow can then be run as usual using `dx run`. For example, if the workflow takes string argument `X`, then: ``` dx run foo -i0.X="hello world" ```

Compilation can be controled with several parameters.

| Option   |  Description |
| ------   | ------------ |
| archive  | Archive older versions of applets.|
| compileMode [All, IR] | Compilation mode - If not specified, the compilation mode is "All" and the compiler will translate WDL or CWL inputs into DNAnexus workflows and tasks.  Use "IR" if you only want to parse CWL or WDL files and convert standard-formatted input files to DNAnexus JSON input format without performing full compilation.|
| defaults \<string\> | JSON file with standard-formatted default values. |
| defaultInstanceType \<string\> | The default instance type to use for "helper" applets that perform runtime evaluation of instance type requirements. This instance type is also used when the '-instanceTypeSelection dynamic' option is set. This value is overriden by any defaults set in the JSON file specified by '-extras'.|
| destination \<string\> | Full platform path (project:/folder) |
| execTree \[json,pretty\] | Print a JSON representation of the workflow. |
| extras \<string\> | JSON file with extra options |
| inputs \<string\> | JSON file with standard-formatted input values. May be specified multiple times. A DNAnexus JSON input file is generated for each standard input file. |
| instanceTypeSelection \[static,dynamic\] | Whether to select instance types at compile time for tasks with runtime requirements that can all be statically evaluated (the default "static" option), or to defer instance type selection in such cases to runtime (the "dynamic" option). Using static instance type selection can save time, but it requires the same set of instances to be accessible during WDL/CWL compilation and during the runtime of the generated applets and workflows. Use the "dynamic" option if you plan on creating global DNAnexus workflows or cloning the generated workflows between DNAnexus organizations with different available instance types. |
| imports \<string\> | Directory to search for imported WDL or CWL files. May be specified multiple times. |
| locked   | Create a locked workflow. When running a locked workflow, input values may only be specified for the top-level workflow. |
| leaveWorkflowsOpen | Leave created workflows open (otherwise they are closed). |
| projectWideReuse | Look for existing applets/workflows in the entire project before generating new ones. The default search scope is the target folder only. |
| reorg    | Reorganize workflow output files. |
| runtimeDebugLevel \[0,1,2\] | How much debug information to write to the job log at runtime. Log the minimum (0), intermediate (1, the default), or all debug information (2, for internal debugging).
| separateOutputs | Store the output files of each call in a separate folder. The default behavior is to put all outputs in the same folder. |
| streamFiles \[all,none,perfile\]| Whether to mount all files with dxfuse (do not use the download agent), to mount no files with dxfuse (only use download agent), or to respect the per-file settings in WDL parameter_meta sections (default). |
| useManifests | Use [manifests](#manifests) files for all workflow and applet inputs and outputs. Implies -locked. |
| waitOnUpload | Whether to wait for each file upload to complete. |

The following common options can also be specified when compiling a workflow.
| Options | Description |
| ------- | ----------- |
| folder \<string> | Platform folder (defaults to '/'). |
| project \<string\> | Platform project (defaults to currently selected project).
| language \<string\> \[ver\] | Which language to use? May be WDL or CWL. You can optionally specify a version. Currently: i. WDL: draft-2, 1.0, and 1.1, and ii. CWL: 1.2 are supported and WDL development is partially supported. The default is to auto-detect the language from the source file.
| quiet | Do not print warnings or informational outputs.
| verbose | Print detailed logging.
| verboseKey \<module\> | Print verbose output only for a specific module. May be specified multiple times.
| logFile \<path\> | File to use for logging output; defaults to stderr.
### Inputs

The `-inputs` option allows specifying a Cromwell JSON [format](https://software.broadinstitute.org/wdl/documentation/inputs.php) inputs file. An equivalent DNAnexus format inputs file is generated from it. For example, workflow [files](https://github.com/dnanexus/dxCompiler/blob/main/test/draft2/files.wdl) has input file

```json
{
  "files.f": "dx://project-aaaa:file-wwww",
  "files.f1": "dx://project-aaaa:file-xxxx",
  "files.f2": "dx://project-aaaa:file-yyyy",
  "files.fruit_list": "dx://project-aaaa:file-zzzz"
}
```

Note that the project ID should always be specified in dx URIs. This will speed up execution time by preventing the need for a more expensive API call to resolve the file.

The command

```console
java -jar dxCompiler-xxx.jar compile test/files.wdl -project project-xxxx -inputs test/files_input.json
```

generates a `test/files_input.dx.json` file that looks like this:

```json
{
  "f": {
    "$dnanexus_link": {
      "id": "file-wwww",
      "project": "project-aaaa"
    }
  },
  "f1": {
    "$dnanexus_link": {
      "id": "file-xxxx",
      "project": "project-aaaa"
    }
  },
  "f2": {
    "$dnanexus_link": {
      "id": "file-yyyy",
      "project": "project-aaaa"
    }
  },
  "fruit_list": {
    "$dnanexus_link": {
      "id": "file-zzzz",
      "project": "project-aaaa"
    }
  }
}
```

The workflow can then be run with the command:

```console
$ dx run files -f test/files_input.dx.json
```

#### CWL Files

In CWL, files have additional fields that necessitate all file inputs being passed using a specially formatted object, rather than a DNAnexus link (i.e. in the applets generated by dxCompiler for a CWL workflow, all CWL `File` inputs are represented using DNAnexus inputs of class `Hash`). The input value is represented in JSON as an object with a special key (`___`) and a value with the following fields:

* `type`: Must be `"File"`.
* `uri`: The `dx://` URI of the file.
* `basename`: The name to use when localizing the file. Optional, defaults to the source file name.
* `contents`: The contents of the file. Optional. If specified, `uri` is ignored and `basename` must be specified. A file is created on the worker having the given basename and contents.
* `checksum`: The file checksum. Optional. If specified, the checkum of the localized file must match or the job will fail with an error.
* `secondaryFiles`: An array of files/directories that must be localized along side the primary file. The is identical in format to a directory listing (see the next section). Secondary files must be listed explicitly (patterns are not allowed).
`format`: An IRI for the file format. See the [CWL specification](https://www.commonwl.org/v1.2/CommandLineTool.html#File). Optional.

Simple example:

```json
{
  "myapp.myfile": "dx://project-xxx:/path/to/file"
}
```

is transformed into 

```json
{
  "myfile": {
    "___": {
      "type": "File",
      "uri": {
        "$dnanexus_link": {
          "id": "file-xxx",
          "project": "project-xxx"
        }
      }
    }
  }
}
```

More complex example:

```json
{
  "myapp.myfile": {
    "class": "File",
    "basename": "foo.txt",
    "contents": "This goes into the file"
  }
}
```

is transformed into:

```json
{
  "myfile": {
    "___": {
      "type": "File",
      "basename": "foo.txt",
      "contents": "This goes into the file"
    }
  }
}
```

which, on the worker, results in a file `foo.txt` being created in the inputs directory with the given contents. This file can be used like any other input file.

#### Directories

Both CWL and the development version of WDL have a `Directory` data type. Although DNAnexus does not treat folders as 
first-class objects, dxCompiler does support `Directory`-typed inputs and outputs, with some caveats.

A folder within a DNAnexus project can be represented in a standard JSON/YAML input file as a URI of the 
following form: `dx://project-xxx:/path/to/folder/` (note that the trailing `/` is required). When this file is passed 
to dxCompiler via the `-inputs` option, it is transformed into DNAnexus input format. 

##### WDL
###### Directory inputs:
In WDL, directories are represented as strings. For example, if the following directory input in a standard JSON input file:

```json
{
  "mytask.dir": "dx://project-xxx:/path/to/folder/"
}
```

is passed to dxCompiler using the `-input` option, it is transformed into the following DNAnexus JSON input file:

```json
{
  "dir": "dx://project-xxx:/path/to/folder/"
}
```

##### Caveats
1. The WDL specification states that a `Directory` input is to be treated as a snapshot of the directory at the time the 
job is executed. To enforce this behavior, at the start of the job the full (recursive) listing of the directory is 
retrieved, and only those files/subfolders are localized to the worker. This means that if a file is added to or removed 
from the directory in the DNAnexus project while the job is running, that change is not reflected in the local copy on 
the worker. However, if the same directory is used in multiple jobs, there is (currently) no way to guarantee that the 
contents are the same between workers. We strongly recommend to enact policies and practices to prevent modification of 
folders that will be used as input to compiled WDL workflows.

2. A second important caveat, which results from the fact that folders are not treated as first-class objects by DNAnexus, is that, if [job reuse](#job-reuse) is enabled, a job that is run with the same folder input as a previous job (and all other inputs the same) will reuse the previous job outputs regardless of whether the contents of the folder have changed. There are two possible solutions:

   * Disable job reuse when running executables with `Directory`-type inputs.
   * Enact policies and practices to prevent modification of folders that will be used as input when job reuse is enabled.

3. Derives from Caveat #2: when Stage-N has Directory outputs which are passed to the inputs (also Directory) of the Stage-N+1, 
the latter stage needs VIEW permissions to the project to scan the output from Stage-N. This is different from how the i/o of the 
File type is passed between the stages.

###### Directory outputs:
WDL 2.0 supports directory outputs. Specify directory outputs in your tasks/workflows by providing an output path as a string. 
DNAnexus platform does not support Directory as a first class objects, therefore Directory outputs are coerced to the Hash type.
Example:

```wdl
task task_with_dir_outs {
  input {}

  command <<<
    mkdir folderoutput
    echo hello > folderoutput/hello.txt
  >>>
  output {
    Directory outdir = 'folderoutput/'
  }
}
```
The `folderoutput/` directory has to be specified **without** `dx://` prefix or project ID. It has to exist on the worker 
(see `command` above). After the job is done, output files will be delocalized to `project-xxx:/OUTPUT_DIRECTORY/folderoutput` 
where `OUTPUT_DIRECTORY` is specified output directory in the `dx run` command.

##### CWL

Directories always have an input class of `Hash`. The value is represented in JSON using a special key (`___`) and a value with the following fields:

* `type`: must be `"Folder"`
* `uri`: the `dx://` URI of the folder
* `basename`: the name to use when localizing the directory (defaults to the folder name if not specified) 
* `listing`: an array of `File` and/or `Folder` objects representing the directory structure. The listing can be nested to any level.

For example, in a standard JSON input file:

```json
{
  "mytask.dir": "dx://project-xxx:/path/to/folder/"
}
```

which, when passed to dxCompiler using the `-input` option, is transformed into the following DNAnexus JSON input file:

```json
{
  "dir": {
    "___": {
      "type": "Folder",
      "uri": "dx://project-xxx:/path/to/folder/"
    }
  }
}
```

CWL provides a mechanism for ensuring reproducibility of jobs that take directory inputs, via the `listing` field. We strongly recommend that CWL users specify the folder listing for each directory input. A job will only be reused if both the folder and the listing are identical. The ordering of the listing is taken into consideration when making the comparison, so the listing must be generated deterministically. The default behavior of dxCompiler when using the `-input` option is to generate input files with full listings for all directories. An example of a folder with a listing is:

```json
{
  "dir": {
    "___": {
      "type": "Folder",
      "uri": "dx://project-xxx:/path/to/folder/",
      "listing": [
        {
          "$dnanexus_link": {
            "id": "file-xxx",
            "project": "project-xxx"
          }
        },
        {
          "___": {
            "type": "Folder",
            "uri": "dx://project-xxx:/path/to/folder/subfolder/",
            "listing": ...
          }
        }
      ]
    }
  }
}
```

In CWL, there is an additional data type available, `Listing`. A listing is similar to a `Folder`, except that it does not have a `uri` and instead must have a `basename` and a `listing`. Importantly, the items in the listing do not need to be from the same source folder. At runtime, a directory of the specified structure is constructed on the worker. If a CWL-style input JSON/YAML file is passed to the `-inputs` option of dxCompiler, a `Directory` input is automatically converted to a `Listing` input if it specifies a `basename` and `listing` but not a `location`or `path`.

For example:

```json
{
  "mywf.mylisting": {
    "class": "Directory",
    "basename": "mydir",
    "listing": [
      {
        "class": "File",
        "location": "dx://project-xxx:/path/to/dir1/file1"
      },
      {
        "class": "File",
        "basename": "file2",
        "contents": "This is my second file"
      },
      {
        "class": "Directory",
        "location": "dx://project-xxx:/path/to/folder1"
      },
      {
        "class": "Directory",
        "location": "dx://project-xxx:/path/to/folder2"
      }
    ]
  }
}
```

is converted into:

```json
{
  "mylisting": {
    "___": {
      "type": "Listing",
      "basename": "mydir",
      "listing": [
        {
          "type": "File",
          "uri": {
            "$dnanexus_link": {
              "id": "file-xxx",
              "project": "project-xxx"
            }
          }
        },
        {
          "type": "File",
          "basename": "file2",
          "contents": "This is my second file"
        },
        {
          "___": {
            "type": "Folder",
            "uri": "dx://project-xxx:/path/to/folder1",
            "listing": ...
          }
        },
        {
          "___": {
            "type": "Folder",
            "uri": "dx://project-xxx:/path/to/folder2",
            "listing": ...
          }
        }
      ]
    }
  }
}
```

which results in the following directory structure being created on the worker:

```
mydir
|_file1
|_file2
|_folder1
| |_...
|_folder2
  |_...
```

### Defaults

The `-defaults` option is similar to `-inputs`. It takes a JSON file with key-value pairs,
and compiles them as defaults into the workflow. If the `files.wdl` worklow is compiled with
`-defaults` instead of `-inputs`

```console
$ java -jar dxCompiler-xxx.jar compile test/files.wdl -project project-xxxx -defaults test/files_input.json
```

It can be run without parameters, for an equivalent execution.

```console
$ dx run files
```

### Extras 

The `-extras` command line option takes an additional JSON file to set or override metadata and runtime attributes of workflows and tasks during compilation. See [Setting DNAnexus-specific attributes in extras.json](#setting-dnanexus-specific-attributes-in-extrasjson) for details on how to write the extras file.

If this is file `extraOptions.json`:

```json
{
    "defaultRuntimeAttributes" : {
      "docker" : "quay.io/encode-dcc/atac-seq-pipeline:v1"
    }
}
```

Then adding it to the compilation command line will add the `atac-seq` docker image to all
tasks by default.

```console
$ java -jar dxCompiler-xxx.jar compile test/files.wdl -project project-xxxx -defaults test/files_input.json -extras extraOptions.json
```

## Describe WDL workflow to obtain execution tree

You can describe a DNAnexus workflow that was compiled by dxCompiler to get an execution tree presentating the workflow. The execution tree will include information on the executables in the workflow (applets and subworkflows). By default, the execution tree is return as JSON. You can supply a `-pretty` flag to return a pretty print. 

To obtain execution tree from a dxCompiler compiled workflow:

1. JSON - [example](./examples/four_levels.exectree.json)

```bash
java -jar dxCompiler-xxx.jar describe <workflow_id> 
```

2. prettyPrint - [example](./examples/four_levels.exectree.pretty.txt)

```bash
java -jar dxCompiler-xxx.jar describe <workflow_id> -pretty 
```

# Task and workflow inputs

WDL assumes that a task declaration can be overriden
by the caller, if it is unassigned, or assigned to a constant.

```wdl
task manipulate {
  Int x
  Int y = 6
  Int z = y + x
  ...
}
```

In the `manipulate` task `x` and `y` are compiled to applet inputs,
where `y` has a default value (6). This allows the applet caller to
override them. Declaration `z` is not considered an input, because it
is assigned to an expression.

In a workflow, similarly to a task, a declaration is considered an
input if it is unassigned or or assigned to a constant. For example,
workflow `foo` has three inputs: `ref_genome`, `min_coverage`, and
`config`. Variable `max_coverage` is not compiled into an input
because it is assigned to an expression. Note that `config` is an
input, even though it is located in the middle of the workflow.

```wdl
workflow foo {
    File ref_genome
    Float min_coverage = 0.8
    Float max_coverage = min_coverage + 0.1

    call GetVersion
    scatter (i in [1,2,3]) {
        call RandCheck { input: ref=ref_genome, seed=i }
    }

    String config = "test"
    ...
}
```

WDL allows leaving required call inputs unassigned, and
specifying them from the input file. For example, workflow `math`
calls task `add`, but does not specify argument `b`. It can then
be specified from the input file as follows: `{ "math.add.b" : 3}`.

```wdl
task add {
    Int a
    Int b
    output {
        Int result = a + b
    }
}

workflow math {
    call add { input: a = 3 }
    output {
       Int result = add.result
    }
}
```

Currently, dxCompiler does not support this feature. However, there is a [suggestion](MissingCallArguments.md) for limited support.

# DNAnexus files as outputs

User can work with DNAnexus files inside non-dockerized WDL tasks using `dx` commands and create tasks or workflows that output DNAnexus files. This can be used to generate multiple outputs without storing them all on the local filesystem. DNAnexus output files need to be in a closed state.

## Example 1: DNAnexus files as outputs
```wdl
task find_fastq_in_folder {
  input {
    String fastq_folder
  }
  command<<<
    project=$DX_PROJECT_CONTEXT_ID
    folder=$project://~{fastq_folder}
    fq=$(dx find data --name "*.fastq.gz" --path ${folder} --norecurse --brief)
    for file_id in $fq
    do
        echo "dx://${file_id}"
    done
  >>>
  output {
    Array[File] fq = read_lines(stdout())
  }
}
```

## Example 2: Uploading local file to DNAnexus

```wdl
task unzip_files {
    input {
        Array[File] zipped_files
    }
    command <<<
        for f in ~{sep=' ' zipped_files}; do
            gunzip "${f}"
            f_unzipped="${f%%.gz}"
            file_id=$(dx upload "${f_unzipped}" --brief --tag unzipped --project ${DX_WORKSPACE_ID} --wait)
            # clean up the unzipped file to save local disk space
            rm "${f_unzipped}"
            echo "dx://${DX_WORKSPACE_ID}:${file_id}"
        done
    >>>
    output {
        Array[File] files_out = read_lines(stdout())
    }
}
```

# Task metadata and runtime

## Task meta and parameter_meta
A WDL task has two sections where metadata can be specified:

* meta: Provides overall metadata about the task
* parameter_meta: Provides metadata for each of the input parameters

Both of these sections allow arbitrary keys and values; unrecognized keys must be ignored by the workflow engine. dxCompiler recognized specific keys in each section that are used when generating the native DNAnexus applets. The purpose of these keys is to provide the same information that can be specified in the [dxapp.json](https://documentation.dnanexus.com/developer/apps/app-metadata) file. 

### meta section

The following keys are recognized by dxCompiler (with slightly different representation than in `dxapp.json`):

* `title`: A short title for the applet. If not specified, the task name is used as the title.
* `summary`: A short description of the applet. If not specified, the first line of the description is used (up to 50 characters or the first period, whichever comes first).
* `description`: A longer description of the applet.
* `developer_notes`: Notes specifically for developers of the task.
* `types`: An array of DNAnexus [types](https://documentation.dnanexus.com/developer/api/data-object-lifecycle/types).
* `tags`: An array of strings that will be added as tags on the generated applet.
* `properties`: A hash of key-value pairs that will be added as properties on the generated applet. Both keys and values must be strings.
* `details`: An object with an arbitrary set of details about the applet. The following keys are specifically recognized and used by the platform:
  * `advancedInputs`
  * `citations`
  * `contactEmail`
  * `contactOrg`
  * `contactUrl`
  * `exampleProject`
  * `repoUrl`
  * `upstreamLicenses`
  * `upstreamUrl`
  * `upstreamVersion`
  * `upstreamProjects`: licenses of the dependency software and packages. The following keys are required to ensure compliance with open-source licenses: `name`, `repoUrl`, `version`, `license`, and `licenseUrl`, while `author` is optional but good to have.
  * `whatsNew`: The task's change log. There are two different formats that are accepted:
    * A (possibly Markdown-formatted) string
    * An array of versions, where each version is a hash with two keys: `version`, a version string, and `changes`, an array of change description strings. This object will be formatted into a Markdown string upon compilation.

The following keys are also recognized but currently unused, as they only apply to DNAnexus Apps (not Applets):

* `categories`: A list of DNAnexus [categories](https://documentation.dnanexus.com/developer/apps/app-metadata#categories-user-browseable-categories)
* `open_source`: Whether the generated app should be open-source
* `version`: The app version

### parameter_meta section

The [WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#parameter-metadata-section) defines a `parameter_meta` section that may contain key value pairs to assoicate metadata with input and output variables. Currently, the following keywords are supported:

- `stream`, indicates whether or not an input file should be streamed. See [here](#Streaming) for more details
- Direct mappings to [inputSpec and outputSpec keywords in dxapp.json](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications):
  - `help` - `description` is also accepted as an alias for `help`; if the parameter definition is a string rather than a hash, the string is used as `help`.
  - `group` - parameter grouping (used in the DNAnexus web UI).
  - `label` - human-readable label for the parameter (used in the DNAnexus web UI).
  - `patterns` - accepted filename patterns (applies to `File`-type parameters only).
  - `choices` - allowed parameter values; currently, this is limited to primitive (`String`, `Int`, `Float`, `Boolean`) and `File` types parameters (and `Array`s of these types), i.e. it is not allowed for `Map` or `Struct` parameters.
  - `suggestions` - suggested parameter values; currently has the same limitations as `choices`.
  - `dx_type` - maps to the `type` field in dxapp.json; can be either a `String` value or a boolean "expression" (see example below). Applies to `File`-type parameters only.
  - `default` - a default value for the parameter. This is ignored if the parameter's default value is defined in the `inputs` section.

Although the WDL spec indicates that the `parameter_meta` section should apply to both input and output variables, currently the `parameter_meta` section is mapped only to the input parameters.

### Streaming

Normally, a file used in a task is downloaded to the instance, and
then used locally (*localized*). If the file only needs to be
examined once in sequential order, then this can be optimized by
streaming instead. The Unix `cat`, `wc`, and `head` commands are of
this nature. To specify that a file is to be streamed, mark it as such
in the `parameter_meta` section. For example:

```wdl
task head {
    File in_file
    Int num_lines

    parameter_meta {
        in_file : "stream"
    }
    command {
        head -n ${num_lines} ${in_file}
    }
    output {
        String result = read_string(stdout())
    }
}
```

File streaming is an optimization, and there are limiting rules to its
correct usage. The file must be accessed only once, in sequential
order, from the beginning. It need not be read to the end. If the task
does not keep this contract, it could fail in unexpected ways.

Some tasks have empty command sections. For example, the `fileSize`
task (below) calculates the size of a file, but does not need to
download it.  In such cases, the input files are downloaded lazily,
only if their data is accessed.

```wdl
task fileSize {
    File in_file

    command {}
    output {
        Float num_bytes = size(in_file)
    }
}
```

## Task Runtime and Hints

A task declaration has a runtime section where [native WDL runtime attributes](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#mandatory-runtime-attributes) can be specified. Besides, several DNAnexus-specific attributes can also be added to this section if you want customize its execution environment and runtime behavior on the DNAnexus platform.

### Instance type
An instance type is chosen by the compiler to satisfy the execution requirments defined by the native runtime attributes. If you wish to choose an instance type from the
DNAnexus [native](https://documentation.dnanexus.com/developer/api/running-analyses/instance-types) list, this can be done by specifying the `dx_instance_type` key
instead. 
For example:

```
runtime {
   dx_instance_type: "mem1_ssd2_x4"
}
```

If you want an instance that has a GPU chipset, set the `gpu` attribute to true. For example:
```
runtime {
   memory: "4 GB"
   cpu : 4
   gpu : true
}
```
CWL provides a different set of runtime attributes used in the [ResourceRequirement](https://www.commonwl.org/v1.2/CommandLineTool.html#ResourceRequirement), and a DNAnexus instance type will be chosen accordingly. Recognized attributes and its mapping to WDL are listed [here](CWL_v1.2.0_to_WDL_v1.md#requirements).

### Additional DNAnexus-specific runtime settings

There are several parameters that also affect the runtime behavior of an applet and can be supplied upon its creation in its [dxapp.json file](https://documentation.dnanexus.com/developer/apps/app-metadata#annotated-example):
* `runSpec.executionPolicy`: Specifies when to try to automatically restart failed jobs, and how many times.  Note that applet's restartOn setting will be overridden by the parent job's executionPolicy setting when launched from a job during workflow execution.
* `runSpec.timeoutPolicy`: Specifies the maximum amount of time the job can run
* `access`: Specifies additional permissions and resources the applet can access
* `ignoreReuse`: Specifies whether to allow the outputs of the applet to be reused

Similarly, these attributes can be specified in the WDL workflow, but their representation there is slightly different than in `dxapp.json`, and all of them should be specified under the `runtime` section of the WDL task. Also note that the runtime section is different than the metadata section when it comes to attribute values - specifically, object values must be prefixed by the `object` keyword, and mapped values must have their keys in quotes.

* `dx_restart`: Either an integer value indicating the number of times to automatically restart regardless of the failure reason, or an object value with the following keys:
  * `max`: Maximum number of restarts
  * `default`: Default number of restarts for any error type
  * `errors`: Mapping of [error types](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification) to number of restarts
* `dx_timeout`: Either a string value that specifies days, hours, and/or minutes in the format "1D6H30M" or an object with at least one of the keys `days`, `hours`, `minutes`.
* `dx_access`: An object with any of the keys:
  * `network`: An array of domains to which the app has access, or "*" for all domains
  * `project`: The maximum level of access the applet has to the host project it was launched in - a string with any [DNAnexus access level](https://documentation.dnanexus.com/getting-started/key-concepts/projects#project-access-levels)
  * `allProjects`: The maximum level of access the applet has to all projects the user has access to - a string with any DNAnexus access level
  * `developer`: Boolean - whether the applet is a developer, i.e. can create new applets
  * `projectCreation`: Boolean - whether the applet can create new projects
* `dx_ignore_reuse`: Boolean - whether to allow the outputs of the applet to be reused
  
### Native DNAnexus executable

You can also specify a native DNAnexus app(let) that will be called as a task at runtime by adding key `dx_app` in the `runtime` section. See [Calling existing app(let)s](#calling-existing-applets) for more details.

## Example tasks with DNAnexus-specific metadata and runtime

### Example 1: grep for pattern in file

```wdl
version 1.0

task cgrep {
    input {
        String pattern
        File in_file
        Int? max_results
    }
    Int actual_max_results = select_first([max_results, 3])

    meta {
        title: "Search in File"
        tags: ["search", "grep"]
        details: {
          whatsNew: [
            { version: "1.1", changes: ["Added max_results", "Switched to WDL v1.0"]},
            { version: "1.0", changes: ["Initial release"]}
          ]
        }
    }

    parameter_meta {
        in_file: {
          help: "The input file to be searched",
          group: "Basic",
          patterns: ["*.txt", "*.tsv"],
          dx_type: { and: [ "fastq", { or: ["Read1", "Read2"] } ] },
          stream: true
        }
        pattern: {
          help: "The pattern to use to search in_file",
          group: "Advanced"
        }
        max_results: {
          help: "Maximum number of results to return",
          choices: [1, 2, 3],
          default: 3
        }
    }

    command <<<
        grep -m~{actual_max_results} '~{pattern}' ~{in_file} | wc -l
        cp ~{in_file} out_file
    >>>

    output {
        Int count = read_int(stdout())
        File out_file = "out_file"
    }

    runtime {
      docker: "ubuntu:latest"
      dx_instance_type: "mem1_ssd1_v2_x8"
      dx_ignore_reuse: true
      dx_restart: object {
          default: 1,
          max: 5,
          errors: object {
              UnresponsiveWorker: 2,
              ExecutionError: 2,
          }
      }
      dx_timeout: "12H30M"
      dx_access: object {
          network: ["*"],
          developer: true
      }
    }
}
```

### Example 2: alignment with BWA-MEM

```wdl
version 1.0

task bwa_mem {
  input {
    String sample_name
    File fastq1_gz
    File fastq2_gz
    File genome_index_tgz
    Int min_seed_length = 19
    String? read_group
    String docker_image = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    Int cpu = 4
    Int memory_gb = 8
    Int? disk_gb
  }

  String genome_index_basename = basename(genome_index_tgz, ".tar.gz")
  String actual_read_group = select_first([
    read_group,
    "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:ILLUMINA"
  ])
  Int actual_disk_gb = select_first([
    disk_gb,
    ceil(2 * (size(genome_index_tgz, "G") + size(fastq1_gz, "G") + size(fastq2_gz, "G")))
  ])

  command <<<
  set -euxo pipefail
  tar xzvf ~{genome_index_tgz}
  /usr/gitc/bwa mem \
    -M \
    -t ~{cpu} \
    -R "~{actual_read_group}" \
    -k ~{min_seed_length} \
    ~{genome_index_basename}.fa \
    ~{fastq1_gz} ~{fastq2_gz} | \
    samtools view -Sb > ~{sample_name}.bam
  >>>

  output {
    File bam = "${sample_name}.bam"
  }

  runtime {
    docker: docker_image
    cpu: "${cpu}"
    memory: "${memory_gb} GB"
    disks: "local-disk ${actual_disk_gb} SSD"
    dx_timeout: "1D"
    dx_restart: object {
      max: 3
    }
  }

  meta {
    title: "BWA-MEM"
    description: "Align paired-end reads using BWA MEM"
    details: {
      upstreamLicenses: "GPLv3"
    } 
  }

  parameter_meta {
    sample_name: {
      label: "Sample Name",
      help: "Name of the sample; used to prefix output files"
    }
    fastq1_gz: {
      label: "FASTQ 1 (gzipped)",
      description: "Gzipped fastq file of first paired-end reads",
      stream: true
    }
    fastq2_gz: {
      label: "FASTQ 2 (gzipped)",
      description: "Gzipped fastq file of second paired-end reads",
      stream: true
    }
    genome_index_tgz: {
      label: "Genome Index (.tgz)",
      description: "Tarball of the reference genome and BWA index",
      stream: true
    }
    min_seed_length: {
      label: "Minimum Seed Length",
      help: "Matches shorter than INT will be missed.",
      group: "Advanced",
      default: 19
    }
    read_group: {
      label: "Read Group",
      help: "(Optional) the read group to add to aligned reads",
      group: "Advanced" 
    }
    docker_image: {
      label: "Docker Image",
      help: "Name of the docker image to use",
      group: "Resources",
      default: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    }
    cpu: {
      label: "CPUs",
      help: "Minimum number of CPUs to use",
      group: "Resources",
      default: 4
    }
    memory_gb: {
      label: "Memory (GB)",
      help: "Minimum amount of memory required",
      group: "Resources",
      default: 8
    }
    disk_gb: {
      label: "Disk Space (GB)",
      help: "Minimum amount of disk space required (in GB); by default this is calculated from the inputs",
      group: "Resources"
    }
  }
}
```

\* Note the comma seperating the members of the objects within meta and paramter_meta

# Calling existing app(let)s

The DNAnexus tools library provides apps for many existing bioinformatics tools, and you may have already developed app(let)s of your own. You may want to use these existing app(let)s rather than rewriting them in WDL. Calling a native app(let) from WDL can be done using a native task wrapper. The dxCompiler `dxni` subcommand is provided to generate native task wrappers automatically. It can generate a wrapper for a specific app(let), all apps, and/or all applets in a specific platform folder. For example, the command:

```console
$ java -jar dxCompiler-xxx.jar dxni -project project-xxxx -folder /A/B/C -output dx_extern.wdl
```

will find native applets in the `/A/B/C` folder, generate tasks for them, and write to local file `dx_extern.wdl`. If an applet has the `dxapp.json` signature:

```
{
  "name": concat,
  "inputSpec": [
    {
      "name": "a",
      "class": "string"
    },
    {
      "name": "b",
      "class": "string"
    }
  ],
  "outputSpec": [
    {
      "name": "result",
      "class": "string"
    }]
}
```

The WDL definition file will be:

```wdl
version 1.0
  
task concat {
  input {
    String a
    String b
  }
  command {}
  output {
    String c = ""
  }
  runtime {
    dx_app: object {
      id: "applet-xxxx",
      type: "applet" 
    }
  }
}
```

The runtime section includes the ID of the app(let) that will be called at runtime.

A WDL workflow can call the `concat` task as follows:

```wdl
import "dx_extern.wdl" as lib

workflow w {
  call lib.concat as concat {
    input: a="double", b="espresso"
  }
  output {
    String result = concat.c
  }
}
```

## Calling apps

To generate WDL calling apps instead of applets, use

```console
$ java -jar dxCompiler.jar dxni -apps only -o my_apps.wdl
```

The compiler will search for all the apps you can call and create WDL tasks for them. The WDL task will look like:

```wdl
version 1.0
  
task concat {
  ...
    
  runtime {
    dx_app: object {
      id: "app-xxxx",
      type: "app"
    }
  }
}
```

You can also use `dx_app.name` rather than `dx_app.id` to specify the app by name, e.g.

```wdl
version 1.0
  
task concat {
  ...
    
  runtime {
    dx_app: object {
      name: "concat_native/1.0.0",
      type: "app"
    }
  }
}
```

## Calling apps and applets using WDL `development`

In version `development` (aka `2.0`), the `runtime` section no longer allows arbitrary keys. Instead, use the hints section:

```wdl
version development
  
task concat {
  ...
    
  hints {
    dnanexus: {
      "app": {
        "name": "concat_native/1.0.0",
        "type": "app"
      }
    }
  }
}
```

## Overriding the native app and applet instance type

By default, when a native app(let) is called it is run using its default instance type. This can be overridden in a native task wrapper just as it can with a regular task:

```wdl
version 1.0
  
task concat {
  ...
    
  runtime {
    dx_app: object {
      name: "concat_native/1.0.0",
      type: "app"
    }
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}
```
### Unsupported overrides
Current version of dxCompiler **does not** support instance override for native app(let)s in the following cases:
1. Specifications of system requirements (RAM/CPU/storage/etc.) is not supported. Use instance names instead e.g. `mem1_ssd1_v2_x4` (see **Global workflow recommendations** for more information)
2. Instance override by string interpolation: e.g. `mem1_ssd1_~{REMAINING_PART_OF_INSTANCE_NAME}` where `REMAINING_PART_OF_INSTANCE_NAME = v2_x4` to construct `mem1_ssd1_v2_x4` when the string is interpolated
3. Dynamic calculation of required RAM/CPU - if the workflow is compiled with `-instanceTypeSelection dynamic`, instance types will not be overridden and default instances of the native executable will be used: 
```wdl
...
Int calc_mem = 30+15
runtime {
  memory: "~{calc_mem} GB"
}
...
```

# Workflow metadata

Similar to tasks, workflows can also have `meta` AND `parameter_meta` sections that contain arbitrary workflow-level metadata. dxCompiler recognizes the following `meta` attributes which are usually specified in [dxworkflow.json](https://documentation.dnanexus.com/developer/workflows/workflow-metadata) and uses them when generating the native DNAnexus workflow:

* `title`: A short title for the workflow. If not specified, the task name is used as the title.
* `summary`: A short description of the workflow. If not specified, the first line of the description is used (up to 50 characters or the first period, whichever comes first).
* `description`: A longer description of the workflow.
* `types`: An array of DNAnexus [types](https://documentation.dnanexus.com/developer/api/data-object-lifecycle/types).
* `tags`: An array of strings that will be added as tags on the generated workflow.
* `properties`: A hash of key-value pairs that will be added as properties on the generated workflow. Both keys and values must be strings.
* `details`: A hash of workflow details. The only key that is specifically recogized is `whatsNew`, and the formatting is handled for workflows the same way as it is for tasks.

The workflow `parameter_meta` section supports the same attributes as the task `parameter_meta` section.

# Setting DNAnexus-specific attributes in extras.json

When writing a DNAnexus applet the user can specify metadata and runtime options through the [dxapp.json](https://documentation.dnanexus.com/developer/apps/app-metadata#annotated-example) file. The dxCompiler equivalent is the _extras_ file, specified with the `-extras` command line option.

The following first-level keys are accepted in the _extras_ file:
* `defaultRuntimeAttributes`: native WDL/CWL runtime attributes
* `defaultTaskDxAttributes` metadata and runtime attributes defaults for tasks
* `perTaskDxAttributes`: metadata and runtime attributes for specific tasks
* `defaultWorkflowDxAttributes`: metadata and runtime attributes defaults for workflows 
* `perWorkflowDxAttributes`: metadata and runtime attributes for specific workflows
* `dockerRegistry`: private registry configuration. See [Private registries](#private-registries])
* `customReorgAttributes`: custom reorganization applet URI and its configuration. See [Adding config file based reorg applet at compilation time](#adding-config-file-based-reorg-applet-at-compilation-time)
* `ignoreReuse`: boolean value indicating whether to disable [job reuse](#job-reuse)
* `delayWorkspaceDestruction`: boolean value indicating whether to delay the destruction of the job's temporary workspace after execution. See [Delay workspace destruction](#delay-workspace-destruction)

If one attribute is specified multiple times, its final value will be retrieved from the following sources and the latter (if exists) will override the former: 
* `defaultTaskDxAttributes`/`defaultWorkflowDxAttributes` in the extras file
* task/workflow metadata/runtime section
* `perTaskDxAttributes`/`perWorkflowDxAttributes` in the extras file

\* Note: the first-level keys in the extras file have been changed to camel case; however, all the old keys (v2.1.0 and earlier) are still recoginzed.

## Default and per-task attributes

The extras file has a `defaultTaskDxAttributes` section where `executionPolicy`, `timeoutPolicy`, and `access` `runSpec` attributes can be set.

```json
{
  "defaultTaskDxAttributes" : {
    "runSpec": {
        "executionPolicy": {
          "restartOn": {
            "*": 3
          }
        },
        "timeoutPolicy": {
          "*": {
            "hours": 12
          }
        },
        "access" : {
          "project": "CONTRIBUTE",
          "allProjects": "VIEW",
          "network": [
            "*"
          ],
          "developer": true
        }
      }
  }
}
```

In order to override the defaults for specific tasks, you can add the `perTaskDxAttributes` section. For example

```json
{  
  "defaultTaskDxAttributes" : {
    "runSpec": {
        "timeoutPolicy": {
          "*": {
            "hours": 12
          }
        }
      }
    },
  "perTaskDxAttributes" : {
    "Add": {
      "runSpec": {
        "timeoutPolicy": {
          "*": {
            "minutes": 30
          }
        }
      }
    },
    "Inc" : {
      "runSpec": {
        "timeoutPolicy": {
          "*": {
            "minutes": 30
          }
        },
        "access" : {
          "project": "UPLOAD"
        }
      }
    }
  }
}
```

will override the default task timeout and `dx_timeout` in the runtime section of tasks `Add` and `Inc`. It will also provide `UPLOAD` instead of `VIEW` project access to `Inc`.

You are also able to specify metadata for tasks in the `defaultTaskDxAttributes` and `perTaskDxAttributes` sections, including adding citation or license information. The full set of recognized attributes is listed in [task meta](#meta-section) with snake_case fields converted to camelCase to agree with [dxapp.json](https://documentation.dnanexus.com/developer/apps/app-metadata) syntax (e.g.  `developer_notes` and `open_source` meta attributes are specified as `developerNotes` and `openSource` in `extras.json`).

For example:

```json
{
  "defaultTaskDxAttributes": {
    "version": "1.0.0"
  },
  "perTaskDxAttributes" : {
    "Add": {
      "runSpec": {
        "timeoutPolicy": {
          "*": {
             "minutes": 30
          }
        }
      },
      "details": {
        "upstreamProjects": [
          {
            "name": "GATK4",
            "repoUrl": "https://github.com/broadinstitute/gatk",
            "version": "GATK-4.0.1.2",
            "license": "BSD-3-Clause",
            "licenseUrl": "https://github.com/broadinstitute/LICENSE.TXT",
            "author": "Broad Institute"
          }
        ]
      }
      "developerNotes": "This note is used to provide additional info to advanced users."
    }
  }
}
```

Note that `details` and `developerNotes` specified in `perTaskDxAttributes` override `details` and `developer_notes` that are set in the task's `meta` section. 

## Default and per-workflow attributes

There are also metadata [attributes](#workflow-metadata) that can be set at the workflow level. You can specify metadata as defaults as `defaultWorkflowDxAttributes`, or set them for each workflow in the `perWorkflowDxAttributes` section.

Note that workflows in WDL do not have a native runtime section; however, you can set the `chunkSize` limit for scatters in the workflow level as a DNAnexus-specific runtime attribute. DNAnexus executes large scatters in "chunks" of no more than 1000 jobs at a time (the default is 500). For some scatters, it may be necessary to increase or decrease the chunk size for efficient execution. You should not need to modify this attribute unless instructed to do so by the DNAnexus support team.

Consider the following workflow:

```wdl
workflow wf1 {
  input {
    Array[Array[File]] samples
    Array[Int] numbers
  }
  
  scatter (sample_files in samples) {
    scatter (file in sample_files) {
      call summarize { input: file = file }
    }
  }
  
  scatter (num in numbers) {
    call add { input: num = num }
  }
  
  output {
    Array[String] summary = mytask.summary
  }
}

task summarize { ... }
task add { ... }
```

If you want the default scatter chunk size for this workflow to be 100, but you want the scatter chunk size for nested scatter (`scatter (file in sample_files) { ...}`) to be 700, then you'd use the following configuration: 

```json
{
  "perWorkflowDxAttributes": {
    "wf1": {
      "scatterDefaults": {
        "chunkSize": 100
      },
      "scatters": {
          "sample_files.file": {
          "chunkSize": 700
        }
      }
    }
  }
}
```

## Job reuse

By default, job results are [reused](https://documentation.dnanexus.com/user/running-apps-and-workflows/job-reuse). This is an optimization whereby when a job is run a second time, the results from the previous execution are returned, skipping job execution entirely. Sometimes, it is desirable to disable this behavior. To do so, add this setting to the top-level of the extras file:

```
{
  "ignoreReuse" : true
}
```
This will be applied to the top-level workflow, sub-workflows, and applets during compilation, and used for all jobs/analyses during execution (which is equivalent to using `--ignore-reuse` flag with `dx run`). 

## Delay workspace destruction

When calling a workflow with `dx run`, jobs and analyses launched by this workflow will have their temporary workspaces to store resources and intermediate outputs. By default, when a job or an analysis has transitioned to a terminal state (done, failed, or terminated), its temporary workspace will be destroyed by the system. 

If you wish to delay the destruction of these temporary workspaces and keep them around longer (around 3 days) after executing the workflow, two things need to be done:

1. Add this setting to the top-level of the extras file:

```
{
  "delayWorkspaceDestruction" : true
}
```
This will guarantee the workspace containers of all jobs that spawned from the parent jobs (e.g. in scatters) during workflow execution to remain intact after the analysis.

2. When running the workflow use `--delay-workspace-destruction` flag:
```
dx run YOUR_WORKFLOW --delay-workspace-destruction
```
This will guarantee the root analysis and its stages will have their workspace destruction delayed.

Taking both steps will ensure all of your workflow containers are not immediately destroyed regardless of whether the analysis succeeds or fails, and will allow you to view the data objects stored in the workflow workspace containers for debugging purposes.

To learn about jobs' workspaces used during execution, please refer to [the official DNAnexus documentation](https://documentation.dnanexus.com/user/running-apps-and-workflows/job-lifecycle#example-execution-tree).

# Handling intermediate workflow outputs

A workflow may create a large number of files, taking up significant disk space, and incurring storage costs. Some of the files are workflow outputs, but many of them may be intermediate results that are not needed once the workflow completes. By default, all outputs are stored in one platform folder. With the `-reorg` flag, the intermediate results are moved into a subfolder named "intermediate". This is achieved by adding a stage to the workflow that reorganizes the output folder, it uses `CONTRIBUTE` access to reach into the parent project, create a subfolder, and move files into it.

## Use your own applet

You may want to use a different applet than the one provided with `-reorg`. To do that, write a native applet, and call it at the end your workflow.

Writing your own applet for reorganization purposes is tricky. If you are not careful, it may misplace or outright delete files. The applet:
1. Requires `CONTRIBUTE` project access, so it can move files and folders around.
2. Has to be idempotent, so that if the instance it runs on crashes, it can safely restart.
3. Has to be careful about inputs that are *also* outputs. Normally, these should not be moved.
4. Should use bulk object operations, so as not to overload the API server. 
   
You must also be aware that the analysis information is updated in the platform's database asynchronously, so the result of calling `dx describe` on the analysis may not be up-to-date. The most reliable method for making sure you have an up-to-date analysis description is to call `dx describe` in a loop (waiting at least 3 seconds between iterations), and exit the loop when the `dependsOn` field returns an array that contains exactly one item - the ID of the reorg job itself. See the [example](CustomReorgAppletExample.md).

## Adding config file based reorg applet at compilation time

In addition to using `-reorg` flag to add the reorg stage, you may also add a custom reorganization applet that takes an optional input by declaring a "customReorgAttributes" object in the JSON file used as parameter with `-extras`

The `customReorgAttributes` object has two properties in extra.json:
* `appUri`: reorg app or applet URI - either an ID (e.g. "app-bwa_mem" or "app-xxx" or "applet-yyy") or a URI of a platform file (e.g. "dx://file-xxx").
* `configFile`: auxiliary configuration file.

The optional input file can be used as a configuration file for the reorganization process.

For example:

```
{
  "customReorgAttributes" : {
    "appUri" : "applet-12345678910",
    "configFile" : "dx://file-xxxxxxxx"
  }
}

# if you do not wish to include an additional config file, 
# you can omit "configFile" or set it to `null`
{
  "customReorgAttributes" : {
    "appUri" : "applet-12345678910",
    "configFile" : null
  }
}
```

The config-file based reorg applet needs to have the following input specs in the dxapp.json:

```json
{
  "inputSpec": [
    {
      "name": "reorg_conf___",
      "label": "Auxiliary config input used for reorganisation.",
      "help": "",
      "class": "file",
      "patterns": ["*"],
      "optional": true
    },
    {
      "name": "reorg_status___",
      "label": "A string from output stage that act as a signal to indicate the workflow has completed.",
      "help": "",
      "class": "string",
      "optional": true
    }
  ]
}
```

When compiling a workflow with a custom-reorg applet declared with `-extras` JSON, a string variable `reorg_status___` with the value of `completed` will be included in the output stage.

The `reorg_status___` is used to act as a dependency to signal that the workflow has completed.

For an example use case of a configuration based custom reorg applet, please refer to [CustomReorgAppletExample.md](CustomReorgAppletExample.md).

# Top-level calls compiled as stages

If a workflow is compiled in unlocked mode, top level calls with no
subexpressions are compiled directly to dx:workflow stages. For
example, in workflow `foo` call `add` is compiled to a dx:stage.
`concat` has a subexpression, and `check` is not a top level call; they
will be compiled to dx:applets.

```wdl
workflow foo {
    String username
    Boolean flag

    call add
    call concat {input: x="hello", y="_" + username }

    if (flag) {
        call check {input:  factor = 1 }
    }
}

task add {
    Int a
    Int b
    command {}
    output { Int result = a + b }
}

task concat {
   String s1
   String s2
   command {}
   output { String result = s1 + s2 }
}

task check {
   Int factor = 3
   ...
}
```

When a call is compiled to a stage, missing arguments are transformed
into stage inputs. The `add` stage will have compulsory integer inputs
`a` and `b`.

For an in depth discussion, please see [Missing Call Arguments](MissingCallArguments.md).

# Manifests

In extreme cases, running compiled workflows can fail due to DNAnexus platform limits on the total size of the input and 
output JSON documents of a job. An example is a task with many inputs/outputs that is called in scatter over a large collection. 
In such a case, you can enable manifest support at compile time with the `-useManifests` option. 
This option causes each generated applet or workflow to accept inputs as an array of manifests, and to produce outputs as a single manifest.

A manifest is a JSON document that contains all the inputs/outputs that would otherwise be passed directly to/from the 
workflow stage. A manifest can be specified in one of two ways: 
1. A `.json` input file (see [Manifest JSON](#manifest-json)) is the recommended way to provide inputs in the manifest format.  
`java -jar dxCompiler.jar -inputs mymanifest.json` will produce `mymanifest.dx.json` that can be passed to `dx run -f mymanifest.dx.json`.
2. A platform `file-xxx` with content described in  [Intermediate manifest file inputs and outputs](#intermediate-manifest-file-inputs-and-outputs) 
section can be used to pass manifest output from a stage of one workflow (including the `output` stage) as input to another workflow.  A 
typical use case for this scenario is when a user wants to pass manifest output file from a stage (including `output` stage) 
directly to a new workflow. Also, this scenario might be useful when debugging individual stages of a failing workflow.

## Manifest JSON

When manifest support is enabled, applet/workflow outputs which are passed from one stage to another (or to the final output
stage) exist in the form of intermediate manifests. Here we describe the format of intermediate manifest for informational purposes only.
There is no need to use them as your workflow inputs, as the JSON manifest above is the recommended format.  
For example, given the following workflow:

```wdl
version 1.1

task t1 {
  input {
    File f
  }
  command <<<
    echo "t1: "  >>out
    cat "~{f}" >>out
  >>>
  output {
     File t1_out = "out"
  }
}

workflow test {
  input {
    String s
    File f
  }
  call t1 {    input:      f = f  }
  output {
    File wf_out = t1.t1_out
  }
}
```

You would write the following manifest:

`mymanifest.json`
```json
{
  "test.input_manifest___": {
    "s": "hello",
    "f": "dx://project-aaa:file-xxx"
  }
}
```

`dxCompiler` also supports input manifests in the format of raw inputs, typically used with Cromwell:
`alternative_manifest.json`
```json
{
  "test.s": "hello",
  "test.f": "dx://project-aaa:file-xxx"
}
```


Compile the workflow `test` from above with the `-inputs mymanifest.json` option. A new file `mymanifest.dx.json` will be 
created with the following content. **NOTE** `mymanifest.dx.json` is created by the compiler - the user does not need to 
create/change it manually.


`mymanifest.dx.json`
```json
{
  "input_manifest___": {
    "encoded": false,
    "types": {
      "f": "File",
      "s": "String"
    },
    "values": {
      "s": "hello",
      "f": "dx://project-aaa:file-xxx"
    }
  }
}
```

The created `mymanifest.dx.json` should be used as an input file when running the workflow:
```commandline
dx run workflow-yyy -f mymanifest.dx.json
```


#### Intermediate manifest file inputs and outputs

When manifest support is enabled, applet/workflow outputs which are passed from one stage to another (or to the final output 
stage) exist in the form of intermediate manifests. Here we describe the format of intermediate manifest for informational purposes only.
There is no need to use them as your workflow inputs, as the JSON manifest above is the recommended format.

Given the above workflow, the manifest output from the `common` stage to the following stages (not shown) would be:

```json
{
  "encoded": false,
  "id": "stage-common",
  "values": {
    "s": "hello",
    "f": "dx://project-aaa:file-xxx"
  }
}
```

The `id` field represents the ID of the stage which created the manifest output. It is optional but will always be 
populated in the output manifests. The manifest may contain additional `types` and `definitions` fields  that are only 
for internal use and can be ignored. The outputs of the workflow are referenced in the `values` field of the output manifest
in the form of a map, where keys are the names of the workflow outputs from the WDL `output` workflow section.

## Analysis outputs

Currently, when a workflow compiled with manifest support is run, the outputs of each job along with the generated 
manifest files are placed directly in the project, in a temporary folder `/.d/<job id>`.  
An example command to extract outputs from the output manifests:  
```bash
ANALYSIS_ID="analysis-yyy"    # your analysis ID here
dx cat $(dx describe --json ${ANALYSIS_ID} | jq -rc '.stages|.[]|select(.id == "stage-outputs")|.execution.output.output_manifest___["$dnanexus_link"]|.project+":"+.id')
```

# Docker

## Setting a default docker image for all tasks

Sometimes, you want to use a default docker image for tasks.
The `-extras` command line flag can help achieve this. It takes a JSON file
as an argument. For example, if `taskAttrs.json` is this file:
```
{
  "defaultRuntimeAttributes" : {
    "docker" : "quay.io/encode-dcc/atac-seq-pipeline:v1"
  }
}
```

Then adding it to the compilation command line will add the `atac-seq` docker image to all tasks by default.

```console
$ java -jar dxCompiler-xxx.jar compile files.wdl -project project-xxxx -defaults files_input.json -extras taskAttrs.json
```

## Private registries

If your images are stored in a private registry, add its information
to the extras file, so that tasks will be able to pull images from it.
For example:

```
{
  "dockerRegistry" : {
      "registry" : "foo.acme.com",
      "username" : "perkins",
      "credentials" : "dx://CornSequencing:/B/creds.txt"
  }
}
```

will `docker login` to `foo.acme.com` with the username of `perkins` and password set to the content of `dx://CornSequencing:/B/creds.txt` prior to fetching docker cointainers.

The credentials are stored in a platform file, so they can be replaced without recompiling. The credentials file must be referenced using a `dx://<project>:<file>` URI, where `<project>` can be a project name or ID, and `<file>` can be a file path or ID. All applets are given the `allProjects: VIEW` permission. This allows them to access the credentials file, even if it is stored on a different project. Care is taken so that the credentials never appear in the applet logs. 

Note that you need to use the full path of the docker image in your WDL. For example, the `myimage:latest` image in the above private registry would be referred to as `foo.acme.com/myimage:latest`.

### AWS ECR registries

Logging into an AWS Elastic Container Registry (ECR) is a bit different than logging into a standard docker registry. Specifically, the AWS command line client is used to dynamically generate a password from an AWS user profile. To handle this use-case, dxCompiler downloads the required AWS `credentials` file, installs the AWS client, and generates the password. See the [AWS documentation](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-profiles.html) for more details and examples.

```
{
  "dockerRegistry": {
    "registry": "<aws_account_id>.dkr.ecr.<region>.amazonaws.com",
    "credentials": "dx://myproj:/aws_credentials",
    "awsRegion": "us-east-1"
  }
}
```

`dx://myproj:/aws_credentials` has AWS credentials:

```
[default]
aws_access_key_id: AKI123ABCDEFT1234567
aws_secret_access_key: ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789/ABC
```

## Storing a docker image as a file

Normally, [docker](https://www.docker.com/) images are public, and
stored in publicly available web sites. This enables reproducibility
across different tools and environments. However, if you have a
docker image that you wish to store on the platform,
you can do `docker save`, followed by uploading the tar ball to platform file `file-xxxx`. Then, specify the docker attribute in the runtime section as
`dx://file-xxxx`. Paths or file ids can be used, for example:
```
runtime {
   docker: "dx://GenomeSequenceProject:/A/B/myOrgTools"
}

runtime {
   docker: "dx://file-xxxx"
}

runtime {
   docker: "dx://project-xxxx:file-yyyy"
}
```

# Proxy configurations

Some organizations place a proxy between internal machines and
external hosts. This is done for security, auditing, and caching
purposes. In this case, the compiler cannot contact the dnanexus API
servers, unless is routes its requests through the proxy. Do
achieve this, set the environment variable `HTTP_PROXY` (or
`HTTPS_PROXY`) to point to the proxy. For example, if you perform the
following on the command line shell:

```bash
$ export HTTP_PROXY = proxy.acme.com:8080
$ java -jar dxCompiler.jar compile ...
```

the compiler will send all requests through the machine `proxy.acme.com` on port `8080`.

If an a proxy with NTLM authentication is used, the following configuration is required:

```
$ export HTTP_PROXY_METHOD=ntlm
$ export HTTP_PROXY_DOMAIN = acme.com
$ export HTTP_PROXY = https://john_smith:welcome1@proxy.acme.com:8080
$ java -jar dxCompiler.jar compile ...
```

# Debugging an applet

## Logging

dxCompiler forwards all of the output (stdout and stderr) from the WDL command to the job log. There is the possibility that excessive logging could cause an out-of-disk-space error. If this occurs, you will need to either use a larger instance type, or reduce the output. To completely ignore output from a command, you can redirect it to `/dev/null`:

```
mycommand > /dev/null 2> /dev/null
```

## Getting applet sources
If you build an applet on the platform with dxCompiler, and want to inspect it, use: ```dx get --omit-resources <applet path>```. This will refrain from downloading the large resource files that go into the applet.

## Getting WDL sources

Compiled workflows and tasks include the original WDL source code in the details field. For example, examine workflow `foo` that was compiled from `foo.wdl`.  The platform object `foo` includes a details field that contains the WDL source, in compressed, uuencoded form. To extract it you can do:

```
dx describe /builds/1.02/applets/hello --json --details | jq '.details | .wdlSourceCode' | sed 's/"//g' | base64 --decode | gunzip
```

# Recompilation

Any significant WDL workflow is compiled into multiple DNAnexus applets and workflows. Naively, any modification to the WDL source would necessitate recompilation of all the constituent objects, which is expensive. To optimize this use case, all generated platform objects are checksumed. If a dx:object has not changed, it is not recompiled, and the existing version can be used. The checksum covers the WDL source code, the DNAnexus runtime specification, and any other attributes. There are two exceptions: the project name, and the folder. This allows moving WDL workflows in the folder hierarchy without recompilation.

# Publishing global workflows

A [global workflow](https://documentation.dnanexus.com/developer/workflows/version-and-publish-workflows#about-workflows-and-global-workflows) is an executable that can be versioned and published to other users. Publishing global workflows may facilitate collaboration across multiple projects, compared with local, project-based workflows.

Publishing a dxCompiler WDL workflow as a global workflow is supported from dxCompiler >= `v2.9.0` and dxpy >= `v0.319.2`. This is done in two steps. First, use `dxCompiler` to compile a workflow from WDL source to a local workflow in a project. Second, use `dx-toolkit` to publish the local workflow as a global workflow. Once the global workflow is published, you can add authorized users.

Example: compiling a WDL workflow for later use as a global workflow.
```
java -jar dxCompiler.jar compile <workflow name>.wdl -instanceTypeSelection dynamic
```

Example: publishing a global workflow from a local workflow. The global workflow's name will match the WDL workflow name. The global workflow's version must be set with `--version`, since a local workflow does not have a `version` property. If `--bill-to` is not specified, your default billing account will be assumed.
```
dx build --globalworkflow --from <project id>:<workflow id> --version <version> --bill-to <user-xxxx | org-yyyy>
dx publish globalworkflow-<workflow name>/<version>
```

Example: [adding and removing authorized users](https://documentation.dnanexus.com/user/helpstrings-of-sdk-command-line-utilities#add-users)
```
dx add users globalworkflow-<workflow name> <user-xxxx | org-yyyy>

dx remove users globalworkflow-<workflow name> <user-xxxx | org-yyyy>
```

Example: [adding and removing tags](https://documentation.dnanexus.com/developer/api/running-analyses/global-workflows#api-method-globalworkflow-xxxx-yyyy-addtags)
```
dx api globalworkflow-<workflow name> addTags '{"tags":["<tag 1>", "<tag 2>"]}'

dx api globalworkflow-<workflow name> removeTags '{"tags":["<tag 1>", "<tag 2>"]}'
```

Example: [adding and removing categories](https://documentation.dnanexus.com/developer/api/running-analyses/global-workflows#api-method-globalworkflow-xxxx-yyyy-addcategories)
```
dx api globalworkflow-<workflow name> addCategories '{"categories":["<category 1>", "<category 2>"]}'

dx api globalworkflow-<workflow name> removeCategories '{"categories":["<category 1>", "<category 2>"]}'
```

Example: [updating title, summary, and/or developer notes](https://documentation.dnanexus.com/developer/api/running-analyses/global-workflows#api-method-globalworkflow-xxxx-yyyy-update)
```
dx api globalworkflow-<workflow name> update '{"title":"<new title>", "summary":"<new summary>", "developerNotes":"<new developer notes>"}'
```

See [Limitations](#limitations) below for more details on which dependencies of the workflow will be automatically included in the global workflow.

## Global workflow recommendations

Avoid storing credentials (passwords, keys, etc.) in the source code of the global workflow, as authorized users will have permission to download (via `dx get`) and view all dxCompiler-generated applets used in the global workflow.

Use simple data types in inputs and outputs of global workflows to make it more intuitive for platform users to provide workflow inputs and examine workflow outputs via CLI and UI.

For better execution stability and to reduce dependence on third-party infrastructure, use Docker images stored on the platform rather than in external registries.

For better portability across projects where the workflow will be run, hard-coding instance types using the key `dx_instance_type` should be avoided for global workflows. You should specify runtime resources using numeric requirements for memory / disk / CPU and compile WDL workflows with the flag `-instanceTypeSelection dynamic`. This option ensures that instance types for jobs will always be selected at runtime, based on the actual instance types available in the runtime project. While this option can result in longer runtimes, it is better for portability because it will never attempt to start a job on an instance type that is not supported.

For informational purposes, include a reference to a git repo commit containing the original workflow source code in the `developerNotes` metadata field (see [example above](#publishing-global-workflows) on how to update developer notes).

Grant appropriate permissions to users authorized to run the global workflow the dependencies of the  global workflow that are not bundled with the global workflow.  These include credentials for external docker registries, DNAnexus apps called within the workflow, and other dependencies discussed in the [Limitations section](#global-workflows-limitations)

## Global workflow limitations

Publishing a dxCompiler-generated workflow as a global workflow is currently only supported for WDL.

The global workflow will currently only support a single region (matching the region in which the original workflow was compiled).

Some dependencies of the original workflow will be automatically included in the global workflow, i.e. they will be cloned into the global workflow's resource container and authorized users of the global workflow will not require additional permissions. These include
- Applets and sub-workflows that were part of the original workflow
- Native applets included in the workflow via `dxni`
- Docker images that are stored as platform files

Some dependencies of the original workflow will not be automatically included in the global workflow, so the user may need additional permissions to access and run the workflow. These include
- Publicly inaccessible DNAnexus apps included in the workflow via `dxni`. Users must have permission to run such apps, which should be granted with `dx add users <app> <user or org>` by apps' developers.
- Platform files referenced in workflow parameters (e.g. default or suggested inputs) or in the workflow body (user needs access to the files)
- Credentials file for a private Docker registry (user needs access to the file)
- Docker images in external registries, or dynamically specified at runtime (these will be pulled at runtime)
- Hard-coded `dx_instance_type` (runtime project needs to support the instance type; using numeric resource requirements, as mentioned under Recommendations, is preferred)

Authorized users will have permission to download (via `dx get`) and view any applets and their data referenced in the global workflow.

<!-- TODO mention URL when the UI supports global workflows -->
Any usage of the above in a workflow (including in its tasks and sub-workflows) will produce a warning in the workflow's `description` metadata field, which can be viewed using:
```
dx describe globalworkflow-<name>/<version> --json | jq -rc '.description | tostring'
```

This also works for a regular workflow:
```
dx describe <project-xxxx>:<workflow-yyyy> --json | jq -rc '.description | tostring'
```
