# Release Notes

## in develop

* Added [DIY Docker image instructions](./scripts/docker_image/). This replaces the Docker image published to Docker Hub, which is deprecated and will be removed after November 1, 2024.

## 2.11.7 2024-09-12

* Update dxda to 0.6.2 and dxfuse to 1.4.0

## 2.11.6 2024-03-26

* Optimize with bulk file describe during scatter collect, also for files inside arrays and structs.

### Dependency updates

#### dxApi [0.13.10](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.10)

## 2.11.5 2024-02-29

* Updated wdlTools and its dependencies
* Updated dxApi
* WDL workflows allow `headJobOnDemand` for the tasks if this spec is passed through `extras` via `runSpec.headJobOnDemand`.

### Dependency updates

#### wdlTools [0.17.17](https://github.com/dnanexus/wdlTools/releases/tag/0.17.17)

* ANTLR4 version bump to 4.13.1

#### cwlScala [0.8.5](https://github.com/dnanexus/cwlScala/releases/tag/0.8.5)

* ANTLR4 version bump to 4.13.1

#### dxApi [0.13.9](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.9)

* adds `headJobOnDemand` attribute to jobNew call

## 2.11.4 2023-07-21
* changes to allow compiling with `treeTurnaroundTimeThreshold` attribute which facilitates platform to send the email 
notifications for the root jobs/analyses with a run time (aka `treeTurnatoundTime`) exceeding the specified threshold. 
The attribute `treeTurnaroundTimeThreshold` can be specified for tasks and workflows via `extras.json`. This feature is 
not exclusive for dxCompiler and more information is available in the platform documentation.

## 2.11.3 2023-06-08
* WDL: eliminates inefficient `file-XXX/describe` API calls and replaces them with a bulk describe (i.e `system/findDataObjects)

### Dependency updates

#### dxApi [0.13.6](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.6)
#### dxFileAccessProtocols [0.5.6](https://github.com/dnanexus/dxScala/releases/tag/protocols-0.5.6)
* changes to facilitate the optimizations of number of `file-xxx/describe` API calls upon (de)localization of input/output files


## 2.11.2 2023-05-30
* fixes awscli dependency
* dxda version update to v0.6.0 to address occasional filesystem errors seen in DNAnexus jobs.
* Fixes for the `-pretty` flag in the `describe` command.
* WDL: fix for the tasks wrapped in a frag with optional inputs. Those inputs get forcibly evaluated to explicit `null` 
and cause inconsistencies with applet inputs. This is fixed and all implicit `null` (`None`) values of optionals are remaining implicit.

## 2.11.0 2023-03-15
* Option `delayWorkspaceDestruction` was deprecated and will be ignored when included in `extras.json`. To preserve the
  workspace, `--delay-workspace-destruction` flag should be included in `dx run workflow-XXX`.
* Minor fixes.

## 2.10.9 2023-03-07

* dxda update to `v0.5.12` which adds a retry for file download. This fixes the "context canceled" error which was
  thrown due to the `dxda` behavior.

### Dependency updates

#### dxApi [0.13.5](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.5)


## 2.10.8 2023-02-23

* WDL: fixes a runtime error when a scatter over an empty array tries to return an empty array. Error occurred 
specifically when a task wrapped in a scatter fragment had an array among the outputs. Currently, evaluation of a task 
over an empty array returns an empty array, and the wrapped task is not executed.
* WDL: fixes a runtime error when a private variable is accessed by index then by key (e.g. `my_var[0].left`).
* Minor: error message now suggests to compile with `-archive` flag instead of the deprecated `-a` flag when an applet 
with duplicate name exists.
* Minor. WDL: removed obsolete warning messages when compiling with `parameter_meta` section in WDL source. 

### Dependency updates

#### dxApi [0.13.4](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.4)

* Fixes handling of `suggestions` field for IO specs when only project ID is specified. Now instead of exception
  a warning is thrown.


## 2.10.7 2022-12-13

* `dxda` update to `v0.5.10` which adds a retry for file download. This fixes the "context canceled" error which was 
thrown due to the `dxda` behavior.
* Added support for new region in Azure London, specifically for OFH projects

### Dependency updates

#### wdlTools [0.17.15](https://github.com/dnanexus/wdlTools/releases/tag/0.17.15)

* Minor fixes and refactoring

## 2.10.6 2022-10-25

* Initial support for `Directory` outputs in the WDL `development` version. Note that the dxCompiler support for the `development` version of WDL is under development and not yet production-ready.
* Minor fixes of type coercion. It includes a fix to the "Value V_String(1.0 GiB) could not be coerced to one of Vector(T_Int)" error when setting `Runtime.memory` as a string in a `development` version task, for example:
```
task {
  ...
  runtime {
    memory: "1 GiB"
  }
}
```

### Dependency updates

#### wdlTools [0.17.14](https://github.com/dnanexus/wdlTools/releases/tag/0.17.14)

* Support for Directory outputs in the WDL `development` version
* Minor fixes of type coercion


## 2.10.5 2022-09-29

* Fix for callable recompilation if the Struct composition was changed

### Dependency updates

#### wdlTools [0.17.13](https://github.com/dnanexus/wdlTools/releases/tag/0.17.13)

* Fix for V_Directory in Struct to Directory input


## 2.10.4 2022-08-23

* CWL: fixed target tool id when overriding docker requirement. File ID of an image provided in `DockerRequirement.dockerLoad` 
in a CWL workflow is now correctly mapped and detected upon override. 
* WDL: manifest mode correctly handles non-fully qualified file IDs (updated dxScala: api)

### Dependency updates

#### dxApi [0.13.3](https://github.com/dnanexus/dxScala/releases/tag/api-0.13.3)

* Handling of non-fully qualified file IDs for bulk search/describe. Now for the files provided without the project ID,
  the `describe` response will be returned only for current workspace/project. If a file was cloned to other projects, they
  will be ignored. Non-fully qualified file IDs are not allowed when searching files in other projects.
* Regression tests for API calls to platform

#### dxCommon [0.11.4](https://github.com/dnanexus/dxScala/releases/tag/common-0.11.4)

* Fix for handling 503 error in dx CLI when API requests are throttled

#### dxFileAccessProtocols [0.5.5](https://github.com/dnanexus/dxScala/releases/tag/protocols-0.5.5)

* Dependency update

#### wdlTools [0.17.12](https://github.com/dnanexus/wdlTools/releases/tag/0.17.12)

* Optional types are preserved when using outputs from the standard library functions without explicit declaration.


## 2.10.3 2022-08-02

* WDL: (WDL >= 1.1) Fix for nested workflows when compiled in the unlocked mode: optional inputs with `None` as default are coerced correctly.  
* CWL: Fix for merging optional source inputs. If an input is a collection (e.g. an array), where some items are the type of `"null"`, it is correctly merged when MultipleInputFeatureRequirement is specified.
* CWL: Fix for making target step argument for cwltool.
* CWL: Added support for "metadata" field for input files in CWL under development (pending inclusion in CWL standard).
* CWL: fix casting cwl types from scatter to scatter.
* CWL: Fix for parameter mapping for output stage of the workflow.

### Dependency updates

#### cwlScala [0.8.4](https://github.com/dnanexus/cwlScala/releases/tag/0.8.4)
* For duplicated calls - prepends a process name to make the name unique, avoiding exception "two different processes with the same name"
* Added support for "metadata" field for input files in CWL under development (pending inclusion in CWL standard).


## 2.10.2 2022-05-17

* WDL: Fix to native app(let)s instance override with system requirements (cpu/memory/disks) in case of the direct calls and executions within fragments. Only `-instanceTypeSelection static` (default) is supported. If compiled with `dynamic` - default instances will be used.
* CWL: Partial implementation of task/fragment applets and workflows reuse. Same as in WDL, the `DocContents` attribute, which is the part of source code that defines each executable, is now used for checksum calculation and comparison. However, since the CWL source code is packed as a nested JSON file (instead of standalone blocks in WDL), fragment applets and workflows will also include the code of all underlying processes in their `DocContents`, and they will not be reused if the wrapped applets/subworkflows have changed.
* CWL: Fix evaluating scatter inputs during workflow execution. It would no longer raise an error if optional inputs of the underlying task are not declared at the scatter level.
* CWL & WDL: Keep the Docker container after a task finishes running, so users can debug Docker related job failures.
* CWL: Fix to the `pattern` handling from `secondaryFiles` at workflow level during compilation

### Dependency updates

#### wdlTools [0.17.11](https://github.com/dnanexus/wdlTools/releases/tag/0.17.11)

* Keep the Docker container after a task finishes running, so users can debug Docker related job failures.
* Runtime object is aware if it was created with default system requirements

#### dxCommon [0.11.3](https://github.com/dnanexus/dxScala/blob/develop/common/RELEASE_NOTES.md#0113-2022-05-11)

* Fix `JsUtils.makeDeterministic` to handle `JsArrays` sorting

#### cwlScala [0.8.3](https://github.com/dnanexus/cwlScala/releases/tag/0.8.3)

* Updated cwljava with fixes to secondaryFiles with pattern field when parsing workflow-level parameters and the helper function utils.Uris.shortname to generate enum symbols without namespaces
* Updated cwljava to 1.0

## 2.10.1 2022-04-18

* WDL: Fragments and blocks reuse applets as well, i.e. they are not rebuilt if the code corresponding to them hasn't been updated in the WDL source file. Previously only tasks were reused. **Breaking**: can break the logic of any App reuse for CWL (even tasks maybe are not reused), because `ApplicationCompiler` and `WorkflowCompiler` now look at the `DocContents` for checksum, and `SourceCode` attribute is now ignored.
* CWL: `NetworkAccess`, `WorkReuse` and `ToolTimeLimit` hints are now supported
* WDL: Update custom reorg applet (used for [custom handling of your workflow outputs](https://github.com/dnanexus/dxCompiler/blob/444036acb16f2555d3cfe5f4c892b9996a8079dc/doc/ExpertOptions.md#adding-config-file-based-reorg-applet-at-compilation-time)) example in documentation.
* WDL: Fix to the order of precedence for the different job reuse settings. The [ignoreReuse](https://github.com/dnanexus/dxCompiler/blob/d6371c3f5087c9de23e671928a741007280c2c33/doc/ExpertOptions.md#setting-dnanexus-specific-attributes-in-extras-file) setting specified in the `extras.json` file should override the [dx_ignore_reuse](https://github.com/dnanexus/dxCompiler/blob/d6371c3f5087c9de23e671928a741007280c2c33/doc/ExpertOptions.md#additional-dnanexus-specific-runtime-settings)  setting specified in the `runtime` section of the WDL file.
* Updates to the documentation, esp. sections about [delay workspace destruction](https://github.com/dnanexus/dxCompiler/blob/444036acb16f2555d3cfe5f4c892b9996a8079dc/doc/ExpertOptions.md#delay-workspace-destruction), [DNAnexus-specific runtime settings](https://github.com/dnanexus/dxCompiler/blob/444036acb16f2555d3cfe5f4c892b9996a8079dc/doc/ExpertOptions.md#additional-dnanexus-specific-runtime-settings), and [outputing DNAnexus files](https://github.com/dnanexus/dxCompiler/blob/444036acb16f2555d3cfe5f4c892b9996a8079dc/doc/ExpertOptions.md#dnanexus-files-as-outputs).
* WDL: Fix to configuring network access to tasks (applets). Now you can disable network access to your tasks with the `runSpec.access` setting in the `extras.json` file:

```
{
  "defaultTaskDxAttributes" : {
    "runSpec": {
        "access" : {
          "network": []
        }
      }
  }
}
```
However, if the task needs network access since it uses a Docker image from an external registry (e.g. DockerHub), the setting will be overwritten and full network access will be given to the task. You can prevent this by storing the Docker image as a file in your DNAnexus project.

Note that the setting above is under `defaultTaskDxAttributes`, which means it will be the default setting for all tasks; follow [these instructions](https://github.com/dnanexus/dxCompiler/blob/444036acb16f2555d3cfe5f4c892b9996a8079dc/doc/ExpertOptions.md#default-and-per-task-attributes) to configure a specific task.

### Dependency updates

#### cwlScala 0.8.1

* Fixes CWL default requirement classnames `NetworkAccess`, `WorkReuse` and `ToolTimeLimit` so the corresponding hints can be recognized by dxCompiler (instead of being defined as `GenericHints` which are not interpreted during compilation).

#### wdlTools 0.17.9

* `TAT.Workflow` has a `source` attribute in analogy to `TAT.Document` to be used for checksum calculation for App/Job reuse
* Fixes evaluation of structs as input parameters of workflow/scatter: hash inputs are coerced from object to struct, and their optional elements are assigned to Null if not specified in job inputs.

## 2.10.0 2022-03-17

* CWL: All tasks and workflows conformance tests now pass 
* CWL: Allow a valid CWL identifier that starts with a number or other characters disallowed by the DNAnexus platform
  * Valid CWL IDs starting with disallowed characters will be padded with an extra "\_\_\_" when generating the full encoded name, which will still be a valid CWL ID (as an URI) when decoded
* CWL: Fix to an error with execution of nested workflow step due to incorrect step identifiers. Previously it prevented running of workflows with the `--single-process` flag.
* WDL: Minor fixes with directory path handling when passing `-separateOutputs` flag
  * When a workflow was [compiled with the -separateOutputs flag](https://github.com/dnanexus/dxCompiler/blob/fa748771fbefe30bb9d311fc891b36696d683eaa/doc/ExpertOptions.md#compiling-workflow), it would result in an error at runtime saying "InvalidInput: Folder path must start with /."
* WDL: Fix an issue with optional identifiers not found when they were a part of an unused expression
* Fix an issue where a `dbcluster` object stored in a project in which a workflow was run was causing a runtime error, if the workflow was compiled with the `-reorg` flag.

### Dependency updates

#### dxApi 0.13.2

* Added `database` and `dbcluster` to the list of DNAnexus data objects that are recognized and can be described

#### dxFileAccessProtocols 0.5.4

* Fix to the issue where `-separateOutputs` option was causing a dx API runtime error due to a missing leading `/` in folder path

#### wdlTools 0.17.8

* Input variables are now evaluated and added to the dependency graph, even when they are optional and not used. They need to be a part of the graph in the case when they are assigned to other (unused) variables for a proper expression evaluation.

## 2.9.1 2022-02-25

* Exception messages of subprocesses are propagated to the exception message of the main process. 
* Fixes a bug where inputs to a WDL block were being ignored if they were Optional
* CWL: run input default is now used when not declared in step
* Minor updates to Developing docs

### Dependency updates

#### dxCommon 0.11.2

* Minor changes to JSON formatting

## 2.9.0 2022-02-08

* Command file is now echoed to stderr rather than stdout
* Fixes a bug where directory support was broken for WDL; directories in WDL are represented as strings
* Switches from using cwltool to cwlpack
* Fix an issue with default values in CWL
* Adds support for publishing global workflows from dxCompiler-generated WDL workflows; see [documentation](https://github.com/dnanexus/dxCompiler/blob/develop/doc/ExpertOptions.md#publishing-global-workflows)
* Fix an error message detecting unsupported CWL version
* Update sbt to 1.6.1
* Update dxda to 0.5.9 and dxfuse to 1.0.0

### Dependency updates

#### cwlScala 0.8.0

* **Breaking**
  * Parser API has changed substantially, with parameters added, removed, and rearranged
  * `Parser.parse` method is now private - use `parseString` or `parseFile` instead
  * Removes all options to modify IDs during parsing
* Handles workflows packed by cwlpack
* Adds `Process.simpleName` method to return simplified process name from ID automatically generated by `cwlpack --add-ids`
* Adds `Identifiable.copySimplifyIds` method to deep-copy objects with simplified IDs
* Adds `CwlEnum.symbolNames` function for getting enum symbols without any namespace prefixes
* Adds `coerce` option to `Evaluator.evaluate` to actually perform type coercion, rather than just checking that the result is coercible to the specified type
* Trims `StringValue` when coercing to primitive types
* Fixes `Evaluator.finalizeInputValue` for compound and optional types
* `Evaluator.finalizeInputValue` loads file contents from remote file source if file does not exist locally
* `CwlType.coerceTo` now returns both the coerced-to type and value
* Added `CwlType.CwlGenericRecord`, which is coercible to either `CwlInputRecord` or `CwlOutputRecord`
* Fixes evaluation of values with multiple possible types

## 2.8.3 2022-01-07

* Update dxCommon and wdlTools - fixes forwarding of stdout/stderr to job log for commands run in docker
* Improves error message when input to call is missing/null

### Dependency updates

#### dxCommon 0.11.1

* Fix `SysUtils.runCommand` forwarding of stderr

#### wdlTools 0.17.7

* Correctly attaches to docker stdout/stderr

## 2.8.2 2022-01-05

* Fixes issue where task without `runtime` section succeeds even when command block results in a failure code
* Forwards command output to the job log rather than buffering it until after command completion
* Updates code to compile with JDK11
* Updates build environment to JDK11, Scala 2.13.7, and SBT 1.5.7

### Dependency updates

#### dxCommon 0.11.0

* Adds `SysUtils.runCommand`, which exposes options for how to handle stdin/stdout/stderr.
* Adds `Paths.BaseEvalPaths.isLocal` attribute to differentiate local from remote paths.
* **Breaking** removes `SysUtils.execScript`. Use `runCommand` with the script path as the argument instead.
* Adds validation of path characters to `FileUtils.getUriScheme`

#### dxFileAccessProtocols 0.5.3

* Handles `ResourceNotFoundException` in `DxFileSource.exists`

#### wdlTools 0.17.6

* Fixes `stderr()` function - previously it was returning the file for stdout

## 2.8.1 2021-12-13

* Excludes apps from `bundledDepends`
* Localizes files declared in WDL task private variables
* Respects runtime definition in native task stub
* Fixes archiving of executables - tags them as "dxCompilerArchived" rather than renaming them
* Fixes error when parsing a field name with multiple `stage-*` prefixes (specifically with stage number >= 10)
* Logs full output of command execution 

### Dependency updates

#### dxCommon 0.10.1

* Adds `PosixPath` class for working with POSIX-style paths
* Changes Paths to use `PosixPath` rather than `java.nio.Path`
* Prettifies truncated log messages

#### dxApi 0.13.0

* Adds `DxApi.addTags` method
* Fixes `DxFindDataObjects` when used with `tags` constraint
* Handles record results in `DxFindDataObjects`
* Adds `systemRequirements` to `DxWorkflowStageDesc`

#### dxFileAccessProtocols 0.5.2

* Uses `PosixPath` rather than `java.nio.Path` for manipulating remote paths

## 2.8.0 2021-11-29

* Indicates whether static instance type selection was used in workflow description annotations and metadata
* Supports cloning workflows between projects, a prerequisite for publishing global workflows
* Fixes dxda manifest downloads for tasks
* Affected by a bug for workflows that include a native platform app via `dxni`. Bug fixed in v2.8.1.

## 2.7.2 2021-11-20

* Uses manifest files for large subjob inputs
* Uses dxda to bulk-download manifest files
* Increases number of retries when downloading single manifest files
* Fixes error when parsing a field name with multiple `stage-*` prefixes
* Allows file-to-string coercion for WDL inputs
* Affected by a bug for workflows that include a native platform app via `dxni`. Bug fixed in v2.8.1.

### Dependency updates

#### dxApi 0.12.0

* Enables `retryLimit` to be set for `DxApi.uploadFile` and `DxApi.downloadFile`

## 2.7.1 2021-11-18

* Fixes issue where jobs fail due to out-of-disk error due to excessive logging
* Adds option to specify native app information in `runtime` section (or `hints` for WDL 2.0)
* Fixes regression where default instance type was overridden when calling a native app
* Affected by a bug for workflows that include a native platform app via `dxni`. Bug fixed in v2.8.1.

### Dependency updates

#### dxCommon 0.9.0

* Adds option to `Logger.trace*` to show beginning and/or end of log when limiting trace length

#### dxApi 0.11.0

* Fixes `DxApi.resolveApp` to handle app name with with version (e.g. `bwa_mem/1.0.0`)
* Adds `version` field to `DxAppDescribe`

#### wdlTools 0.17.4

* Fixes `eval.Meta.get` to respect override values

## 2.7.0 2021-11-10

* Handles empty scatters in WDL workflows
* Reserved parameters are now placed in the "Reserved for dxCompiler" parameter group (only affects the display of the app/workflow in the UI)
* Enables other metadata (title, description, version, etc.) to be set via extras.json
* When using manifests, passes any expression values from the helper applet to the called applet workflow so they are added to the output manifest
* When an applet calls executables, adds executables to `bundledDepends` to support cloning of workflows
* Executor creates hard- rather than soft-links in the input directory, so that linked files are accessible from within containers
* Fixes bug where default input values were not overridden for task inside subworkflow
* Adds support for specifying native app(let) in `runtime` section
* Fixes some type conversion bugs related to CWL `Any` type
* Affected by a bug for workflows that include a native platform app via `dxni`. Bug fixed in v2.8.1.
* Logs entire contents of WDL command at runtime

*Warning*: we discovered a regression in this release that may cause tasks to fail with out-of-disk errors due to excessive logging. Please update to 2.7.1 or later.

### Dependency updates

#### dxApi 0.10.1

* Fixes `resolveProject` to handle `container-` objects
* Improves error message when API call fails due to connection error

#### cwlScala 0.7.1

* Handles error when listing folder during path value finalization

#### wdlTools 0.17.2

* Fixes infinite loop when calling `wdlTools.eval.Runtime.contains` with "docker" or "container"

## 2.6.0 2021-10-11

* Implements `Directory` support for CWL and for WDL development/2.0
* Implements CWL `secondaryFiles` support
* Implements CWL workflow support
* **Breaking Change**: the inputs folder is now modified to be read-only for non-streaming files (this has always been the case for streaming files). This means that tasks that create files in the inputs folder will no longer work. For example, a task that creates an index file in place for an input BAM file should be changed from:
    ```wdl
    command <<<
    samtools index ~{mybam}
    >>>
    ```
  to
    ```wdl
    command <<<
    mkdir bams
    ln -s ~{mybam} bams/~{basename(mybam)}
    samtools index bams/~{basename(mybam)}
    >>>
    ```

### Dependency updates

#### wdlTools 0.17.1
* Fixes parsing of placeholder options in draft-2 and 1.0 such that `default` and `sep` are no longer treated as reserved words

#### cwlScala 0.7.0
* **Breaking Change**: `Sink.linkMerge` is now `Option`al
* Introduces `ParserResult` class, which is returned from all `Parser.parse*` methods
* For packed workflows, parses out `$schemas` and `$namespaces`
* Updates to latest cwljava, which fixes several parsing errors
* Fixes `CwlType.flatten` to correctly handle duplicate types
* Treats scatter sources as identifiers
* Automatically renames the `main` process if its name collides with another process
* Fixes evalution of compound parameter references
* Set `CommandLineTool.successCodes` to `Set(0)` if not specified
* Fixes deserialization of optional fields
* Uses the source file name as the process name when processing a `$graph` where the top-level element ID is 'main'
* Fixes evaluation for Directory-type values with listings
* Fixes parsing of LoadListingEnum values
* Adds option to `Parser.parseFile` and `Parser.parseString` to specify that the CWL file is in "packed" form
* Updates to latest cwljava, which fixes parsing of anonymous schemas in packed documents
* Correctly handles identifiers with namespaces from imported documents
* Fixes error when trying to finalize a File value without location or path
* Updates dxCommon to 0.2.15-SNAPSHOT
* Uses `FileSource.listing` to determine directory listing during finalization
* Updates to dxCommon 0.7.0
* *Breaking change*: schema types now have `id: Option[Identifier]` rather than `name: Option[String]`
* Parser can now handle `$graph` style CWL documents
* Adds dependency on `dxCommon` library  
* Improves finalization of file values
* Fixes coercion of StringValue to CwlEnum
* other bugfixes
* Adds `EvaluatorContext.createInputs` to create an `EvaluatorContext` from input values
* Performs "finalization" of input values (setting of missing attributes on File and Directory values) when using `EvaluatorContext.createInputs` or `EvaluatorContext.createStaticInputs`
* Parser bugfixes
* Incorporate `cwljava/39`, which fixes workflow parsing issues
* Allow duplicate Requirements/Hints and specify priority rules
* *Breaking change*: added new `CwlMulti` type and removed all uses of `Vector[CwlType]`

## 2.5.0 2021-09-14

* Docker image dependencies that are DNAnexus platform files are included in applet's bundledDepends
* Adds `-instanceTypeSelection` compiler option to allow disabling compile-type instance type selection
* Adds `-defaultInstanceType` option
* Adds support for new London region, `aws:eu-west-2-g`
* Fixes `-projectWideReuse`
* Fixes evaluation of `outputs` when one declaration depends on another

### Dependency updates

dxCommon 0.8.0
* Adds `getTargetDir` methods to `LocalizationDisambiguator`
* Fixes use of `localizationDir` together with `force` in `SafeLocalizationDisambiguator`
* Adds `FileSource.resolveDirectory` as a separate method from `resolve`

dxApi 0.10.0
* Removes price-based selection of instance types in favor of rank-based selection
* Fixes parsing of non-file-type default values that are reference-type links
* Fixes parsing of parameter defaults/suggestions/choices that are of type `Hash`

dxFileAccessProtocols 0.5.0
* Implements `resolveDirectory` method

wdlTools 0.17.0
* **Breaking Change** `Eval.applyMap` is changed to `Eval.applyAll` and takes a `Vector` rather than `Map` argument. This is done to ensure the expressions are evaluated in order in case there are dependencies between them.
* Fixes parsing of `runtime.returnCodes`

## 2.4.10 2021-08-16

* Fixes issue with using `File` declarations in scatter/conditional blocks
* Updates dxda to 0.5.7
* Handles empty scatters

## 2.4.9 2021-08-09

* Files that are generated by calling functions in worklfow expressions are now uploaded to the platform
* `dxni` now makes optional any input parameter for which the native app(let) has a default value 
* Jobs now fail on out-of-disk errors
* Updated dxda and dxfuse to latest versions
* Updated dxCommon and dxApi dependencies

## 2.4.8 2021-07-13

* Fixes issue where an optional variable inside a conditional could have an inferred type of `T??`, which is illegal and results in a runtime error
* Specifies operating system when selecting optimal instance type
* Fixes issue with references between output variables
* Uses an output stage if there is an output parameter with a literal value (DNAnexus output parameters do not support default values)
* Fixes issue where compiling with `-execTree [pretty|json]` did not print the workflow tree
* Adds support for specifying app metadata in the hints section in WDL development version
* All compile-time calls to `findDataObjects` now search the context project for any file with no project specified
* All runtime `dxAPI.describeFilesBulk` calls, which invoke the platform's `system.findDataObjects`, now search in the job/analysis workspace first before searching the source project(s). 
  * The call is used to search and  describe input files for the analysis and to search and describe output files during the outputs reorganization of the analysis (in its default implementation, not the customized one). If the files cannot be found in the workspace container they are looked for in the project that qualifies file ID (e.g. in project_A if file_B's path is project_A:file_B)
* Better handles insufficient permissions when requesting instance type price list

### Dependency updates 

#### dxApi [0.6.0](https://github.com/dnanexus/dxScala/blob/release-2021.07.12/api/RELEASE_NOTES.md#dxapi)

* Adds option to `DxApi.describeFilesBulk` to search first in the workspace container
* `DxApi.resolveDataObject` now searches in the current workspace and/or project if the project is not specified explicitly. The call is used to find one data object at a time, it's not used in bulk resolution.
* Refactors` DxFindDataObjects` to use separate `DxFindDataObjectsConstraints` class for specifying constraints
* Uses the currently select project ID as the workspace ID when not running in a job
* Better handles insufficient permissions when requesting instance type price list

## 2.4.7 2021-06-09

* Fixes issue with using struct types in workflow outputs
* Fixes issue with array with optional item type

## 2.4.6 2021-05-27

* Fixes regression in WDL code generator where long expressions within placeholders are line-wrapped incorrectly

## 2.4.5 2021-05-25

* An applet that contains multiple scatter or conditional blocks will now have a name that is the concatenation of all the block names 
* Fixes multiple issues with WDL code generator
* Fixes issue with referencing struct fields/call outputs in declarations within nested blocks

## 2.4.4 2021-05-10

* Escapes WDL strings in generated code
* Fixes issues with using expressions in placeholder option values
* Fixes error when evaluating array element access for an optional value
* Fixes localization of files in identically named folders in different projects
* Fixes localization of files with the same name in the same folder
* Encodes/decodes `dx://` URIs to handle project/file names with spaces

## 2.4.3 2021-04-23

* Fixes an issue where a file input from an different project than where the workflow is compiled is localized to an invalid path
* File downloads (including Docker images) no longer fail when retried after a previous failure
* Adds the `-waitOnUpload` compiler option, which causes all file uploads to block until they complete
* Fixes an issue where tasks with outputs that are collections of files (e.g. `Array[File]`) are compiled with an incorrect default input value
* Fixes an issue where using a field of a struct as a call input causes a runtime error, e.g.
    ```wdl
    struct MyStruct {
      String s
    }
    workflow wf {
      input {
        MyStruct my
      }
      call mytask { input: s = my.s }
    }
    ```

## 2.4.2 2021-04-20

* `-imports` now works correctly - `import` statements can reference files that are relative to an import directory rather than the main document
* Coercion from `Object` to `Map` is now allowed - this enables the result of `read_json` to be assigned to a variable of type `Map`
* Strings with escape sequences are now processed correctly - to preserve escape sequences in a string, they need to be "double-escaped"
  - For example, to pass a regular expression containing a tab character to the second argument of the `sub` function:
      ```wdl
      String s1 = "hello\tBob"
      String s2 = sub(s1, "\\t", " ")
      ```
  - Another common example is passing a read group to `bwa mem`. Whether the string is defined in WDL, in a JSON input file, or on the command line, it needs to be "double-escaped":
      ```wdl
      String rg = "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:ILLUMINA"
      
      command <<<
      bwa -R "~{rg}" ...
      >>>
      ```
* Fixes an issue where indentation is stripped out of some commands
* The default reorg applet now polls the analysis until it is fully updated before proceeding. For custom reorg applets, see the updated code in the [example](doc/CustomReorgAppletExample.md).
* Fixes an issue with nested scatters that reference private variables from outer scopes, e.g.
    ```wdl
    workflow wf {
      input {
        Array[String] samples
        Array[File] files
      } 
      scatter (i in range(length(samples))) {
        scatter (j in range(length(files))) {
          call mytask { 
            input: sample = samples[i], file = files[j]
          }
        }  
      }
    }  
    ```

## 2.4.1 2021-03-25

* Fixes AWS ECR issues: bundles AWS CLI with executor rather than installing at runtime

## 2.4.0 2021-03-18

* Adds `-useManifests` option to generate applets and workflows whose inputs and outputs are manifest files
* Fixes issue with using both streaming and non-streaming file inputs in the same task
* Fixes issue with scatter as the first element of a workflow
* Updates to wdlTools 0.12.7, which provides compatibility for some non-complaint syntax allowed by Cromwell
* Fixes `-separateOutputs` for scatter jobs

## 2.3.1 2021-03-03

* Fixes issue with call arguments that access fields of private declarations
* Fixes issue with passing null to optional call parameters
* Fixes common applet naming issue with complex nested workflows
* **Experimental**: * adds `-useManifests` option to generate applets and workflows whose inputs and outputs are manifest files

## 2.3.0 2021-02-23

* Adds `-separateOutputs` option to store output from each call in a separate folder
* Fixes issue with referencing optional variables in the command block
* Ignores default values from native app stubs that can cause errors during compilation

## 2.2.1 2021-02-17

* Fix: native tasks now use instance type from native app/let, not wrapper task

## 2.2.0 2021-02-12

* DxNI handles native app(let)s with optional non-file object inputs
* Sevaral additions/changes to extras.json:  
  * Adds support for Amazon ECR repositories
  * Adds support for configuring the chunk size for scatters (how many scatter jobs run concurrently)
  * All top-level keys are now camel-case (old-style names are still recognized)

## 2.1.1 2021-02-09

* Adds support for WDL v1.1
* Fixes errors due to expression evaluator not handling values wrapped in V\_Optional
* Fixes errors related to compound references
* Updated wdlTools 0.12.3. This update enforces uniqueness of variable names within the same scope
* Avoids unnecessarily generating sub-workflow output applets
* Fixes issue with task command blocks that begin with a placeholder
* Upgraded dxfuse to 0.24.0

## 2.1.0 2021-02-04

* Update CWL parser to fix compilation of WDL tools with imports
* Add CWL tool compilance tests to integration test suite
* Fixes a bug where structs were out of order in generated code. This bug was due to an unnecessary conversion of a Vector to Map in CodeGenerator, which disordered the items.
* Fixes an issue with resolving nested field references. There were two issues here: i. Unnecessarily including outputs that are direct pass-throughs of inputs in the closure when determining the inputs to an Output stage. ii. A bug in the resolution of deeply nested field references in wdlTools.eval.Eval that is fixed in 0.11.18.
* Added an implementation of the manifest format for WDL.
* Update wdlTools to 0.12.1, dxCommon to 0.2.5, dxApi to 0.1.8, cwlScala to 0.3.4

## 2.0.1 2021-01-22

* Implement translator and executor for CWL (tools only)
* Fix map keys and values being out of order
* Fix to EvalException: identifier not found
* Fix to passing non-optional empty arrays
* Update and fix documentation and WDL examples
* Upgrade [dxfuse](https://github.com/dnanexus/dxfuse/releases) to 0.23.3, which now mounts in a read-only mode by default

## 2.0.0 Release Notes Summary

A summary of changes in the 2.0.0 version of dxCompiler comparing to the dxWDL version 1.47.2 (after which the first 2.0.0-rc version was introduced):

User-facing changes:
* Changed the name and all references from dxWDL to dxCompiler, dxCompiler-\*.jar is now used to compile an application or a workflow
* Replaced Cromwell WOM (parser, type checker and evaluator) with [wdlTools](https://github.com/dnanexus/wdlTools), a library maintained by DNAnexus
* Added support for per-workflow and per-scatter chunk size settings to extras.json
* Added `streamFiles` compile option
* Updated the mechanism of comparing instance types to select the cheapest one for execution
* Improved file resolving and caching
* Optimized bulk description of files by replacing `system/describeDataObjects` with `system/findDataObjects` API call and scoping file search to projects
* Increased file name disambiguation directory size limit to 5000
* Increased number of retries for the DNAnexus API requests to 10

Development and codebase changes:
* Reorganized dxCompiler code into subprojects
  * core: shared code between front-end and back-end
  * compiler: the front-end
  * executorCommon: shared code between executors
  * executorWdl: the WDL executor
  * executorCwl: the CWL executor
* Separated the compiler from executor packages
* Moved API code to a separate dxApi library
* Moved DxFileAccessProtocol to a separate protocols library
* Moved parts of the code to a separate dxCommon library
* Moved WDL type serialization code to wdlTools
* Added archive format and implemented UDFs (user-defined functions) for archiving and unarchiving inputs
* Changed the way the applet\_resources folder is used: common binaries are stored in sub-folders by version, each executor is built in its own sub-folder
* Replaced Travis with Github Actions for unit and integration testing

## 2.0.0 2020-12-17

* Reorganizes dxCompiler code into subprojects, splits codebase into 5 subprojects:
  * core: shared code between front-end and back-end
  * compiler: the front-end
  * executorCommon: shared code between executors
  * executorWdl: the WDL executor
  * executorCwl: the CWL executor
* A few classes were moved into different packages to clear up invalid dependencies (e.g. code in core depending on a class in compiler)
* Changes (almost) all references from dxWDL to dxCompiler
* Updates build and test scripts to handle the new layout
* Changes the way the applet\_resources folder is used: common binaries are stored in sub-folders by version, each executor is built in its own sub-folder
* Updates to latest wdlTools (0.11.8)
* Adds archive format
* Implements udfs for archive and unarchive
* Updates to dxApi 0.1.6 - fixes parsing of DxFindDataObjects workflow results
* Fixes errors encountered for repeated compilation
* Skips epilog if there are no applet outputs
* Fixes conversion of Map schema to WDL Map type - unwrap key and value array types
* Fixes comment parsing bug in WDL v2
 
## 2.0.0-rc6 2020-11-20

- Upgrades wdlTools to 0.11.0
- Fixes the issue with type-checking struct-typed objects
- Adds hooks for user-defined functions, which is needed for the new archive/unarchive operations that will be added
- Uses file name when creating names for scatter jobs
- Uses dnanexus-executable.json to get type information about applet input/output parameters, and uses it to type-check the actual inputs/outputs. Ignores hash types, which might be JSON objects or schemas
- Adds `streamFiles` compile option
- Moves API code to dxApi library
- Moves DxFileAccessProtocol to separate protocols library
- Moves WDL type serialization code to wdlTools
- Makes array type-checking less strict (to coerce between empty and non-empty array types)

## 2.0.0-rc5 2020-11-13

- Upgrades wdlTools to 0.10.5
- Adds dxCommon dependency
- Sets either a default or null value for missing optional block inputs, and throws an exception if there are any missing required block inputs
- Addresses the situation where a call input references a field of the scatter variable
- Updates the mechanism of comparing instance types to select the cheapest one for execution
- Fixes to the "x appears with two different callable definitions" compilation error
- Increases disambiguation dir limit to 5000
- Simplifies job names

## 2.0.0-rc4 2020-10-15

- Upgrade wdlTools to 0.6.1
- Upgrade dxda to 0.5.4
- Increase number of API retries to 10
- Fix to the "x appears with two different callable definitions" compilation error
- Add missing required fields to findXXX/describe
- Fix name regexp in DxFindDataObjects
- Additional fixes and improvements

## 2.0.0-rc3 2020-09-16

- Major code reorganization to separate the compiler from executor
- Upgrade of dxda to v0.5.4 and dxfuse - to v0.22.4

## 2.0.0-rc2 2020-07-30

- Mostly internal changes and code re-organization, in preparation for adding CWL support

## 2.0.0-rc 2020-06-25

- Replaced WOM with `wdlTools`
- TODO:
    - Publish wdlTools to MavenCentral and update `build.sbt`
    - Map `parameter_meta` output parameters, remove note in ExpertOptions
    - Update `Internals.md`
- Replaced Travis with Github Actions for unit testing
- Optimized bulk description of files by replacing `system/describeDataObjects` with `system/findDataObjects` API call and scoping file search to projects

## 1.47.2 2020-06-05
- Upgrade dx-download-agent (includes a fix to the early database close issue) 
- Fix to describing billTo of a project

## 1.47.1 2020-05-26
- Log dxda and dxfuse version in applet execution

## 1.47 2020-05-14
- Improvements to `exectree` option
- Upgrade dxfuse to v0.22.2
- Bug fix in dx-download-agent (https://github.com/dnanexus/dxda/issues/34)

## 1.46.4 2020-04-08
- Limit scatters to 500 elements. Running more than that risks causing platform problems.
- Uprade dxfuse to v0.22.1

## 1.46.3 2020-03-27
- fixed bug when describing a live (non-archived) hidden file.
- Uprade dxfuse to v0.21

## 1.46.2 2020-03-18
- Do not use any of the test instances, or any instance with less than 2 CPUs and 3GiB or RAM for auxiliarly WDL jobs.

## 1.46.1 2020-03-17
- Do not use AWS nano instances for auxiliarly WDL jobs. They are not strong enough for the task.

## 1.46 2020-03-12
- Recognize DNAnexus-specific keys in task and workflow metadata
- Recognize workflow parameter metadata
- Optionally load task and workflow descriptions from README files
- Updated dx-download-agent that reduces memory consumption. This is noticible on small instances with a limited amount of memory.

## 1.45 2020-03-03
- Upgrade packages to:
  - dxfuse v20
  - dx-download-agent with support for symbolic links
  - Cromwell v49
- Retry `docker login` if it fails
- Allow DxNI to work with an [app](https://github.com/dnanexus/dxWDL/issues/364)
- DNAx symbolic links can now be used as input files, the dx-download-agent is able to download them. You will need to specifically allow network access to applets that use symbolic links, otherwise they will won't be able to reach external URLs. For example, the `extras.json` file below sets the timeout policy to 8 hours and allows network access.

```
{
  "default_task_dx_attributes" : {
    "runSpec": {
      "timeoutPolicy": {
        "*": {
          "hours": 8
        }
      },
      "access" : {
        "network": [
          "*"
        ]
      }
    }
  }
}
```

## 1.44 2020-02-21
- Added support for additional parameter metadata:
  - group
  - label
  - choices
  - suggestions
  - dx_type
- Fixed bug in project-wide-reuse option.

## 1.43 2020-02-11
- Providing a tree structure representing a compiled workflow via `--execTree pretty`
- For the json structure use `--execTree json`
- Added support for the parameter_meta: patterns to take an object: [docs/ExpertOptions](doc/ExpertOptions.md#parameter_meta-section)
- Support the `ignoreReuse` and `delayWorkspaceDestruction` in the extras file. More details are in the [expert options](doc/ExpertOptions.md#job-reuse).


## 1.42 2020-01-23
- Providing a JSON structure representing a compiled workflow. This can be done with the command line flag `--execTree`. For example:

```
java -jar dxWDL-v1.42.jar compile CODE.wdl --project project-xxxx --execTree
```

- Bug fix for case where the number of executions is large and requires multiple queries.

## 1.41 2020-01-21
- Added support for `patterns` and `help` in the `parameter_meta` section of a WDL task. For more information, see [docs/ExpertOptions](doc/ExpertOptions.md#parameter_meta-section)
- Upgrade to Cromwell v48
- Upgrade to dxfuse v0.17
- Added support for a `path` command line argument to `dxni`.

```
java -jar dxWDL-v1.41.jar dxni --path /MY_APPLETS/assemble --project project-xxxx --language wdl_v1.0 --output headers.wdl
```

- Using `scalafmt` to normalize code indentation.
- Check if a file is archived prior to downloading or streaming it. Such a file cannot be read, and will cause a
"403 forbidden" http error code.


## 1.40 2019-12-19
- Replaced the dxjava package with scala code. The dnanexus calls now go through the low-level DXAPI java
module.


## 1.37.1 2019-12-18
- Bug fix for calling a task inside a scatter with an optional that is not provided.

## 1.37 2019-12-16

- Fix issue with invalid WOM when compiling workflows containing sub-workflow with expression in output block while using custom reorg applet.
- Warning message for custom reorg applet will only show when `--verbosity` is set.
- Minor changes.

## 1.36.1 2019-11-22

- Upgraded to dxfuse version [0.13](https://github.com/dnanexus/dxfuse/releases/tag/v0.13)
- Making dxfuse startup script more robust

## 1.36 2019-11-18
Added a mechanism for a custom reorganization applet, it can be used instead of the built in --reorg option.
You can use it to reorganize workflow file results after it completes.

Please refer to [docs/ExpertOptions.md](doc/ExpertOptions.md#adding-config-file-based-reorg-applet-at-compilation-time)

## 1.35.1
- Object form of streaming syntax. This allows several annotations for an input/output
parameter. In addition to the previously supported:
```
parameter_meta {
  foo: stream
}
```

You can also write:
```
parameter_meta {
  foo: {
    stream: true
  }
}
```
- New version of dxfuse (v0.12)
- The checksum of an applet/workflow does not include the dxWDL version. This means that upgrading to
a new dxWDL version does not require recompiling everything.

## 1.35 2019-10-24
- Support for GPU instances.

If you want an instance that has a GPU chipset, set the `gpu` attribute to true. For example:
```
runtime {
   memory: "4 GB"
   cpu : 4
   gpu : true
}
```


## 1.34 2019-10-17
- Bug fix for handling of structs when creating task headers

## 1.33 2019-10-15
- Upgrade to Cromwell 47
- Protect 300MiB of memory for dxfuse, if it is running. We don't want it killed by the OOM if the user
processes use too much memory. This works only when using docker images.
- Fixed bug with comparision of tasks.

## 1.32 2019-10-04
- Upgrade to Cromwell 46.1
- Retrying docker pull at runtime.
- Improved release script. Copies to geographically distributed regions with an app.
- Fixed bug [#313](https://github.com/dnanexus/dxWDL/issues/313).

## 1.31 2019-09-30

- Renaming dxfs2 to [dxfuse](https://github.com/dnanexus/dxfuse). This
  is the official name for the DNAx FUSE filesystem.
- [Prefer](https://github.com/dnanexus/dxWDL/issues/309) v2 instances over v1 instances.
- Improving marshalling of WDL values to JSON.
- Allow a WDL input that is a map where the key is a file.

## 1.30
- Using dxfs2 to stream files. This replaces `dx cat`, which was the previous solution.

## 1.22
- Upgrade to Cromwell version 46
- Removed a runtime assert that was too strict. It was checking that the type of a WDL value *V* had a static WDL type *T*. However, the real question was whether *V* could be casted type *T*.

## 1.21.1
- Writing out a better description when raising a runtime assertion because WOM types don't match.

## 1.21
- Documented syntax limitation for task/workflow [bracket placement](README.md#Strict-syntax).
- Upgraded to sbt 1.3.0
- Improvements to the download agent

## 1.20b
- The default timeout limit for tasks is 48 hours. It can be overriden by setting a different timeout in the [extras file](./doc/ExpertOptions.md#Setting-dnanexus-specific-attributes-for-tasks).
- Fixing a bug where the database is locked, in the [dx-download-agent](https://github.com/dnanexus/dxda/).

## 1.20a
- An experimental version of dx-download-agent
- Upgrade to Cromwell 45.1

## 1.20
- Experimental version of the download-agent. Trying to reduce download failures at the beginning of a job. For internal use only.

## 1.19
- Bug fix for DxNI error
- Limit the size of the name of jobs in a scatter
- Upgrade to Cromwell v45
- Nested scatters are supported

## 1.18
- Correctly identify WDL calls with no arguments
- Dealing with the case where a WDL file imports a file, that imports another file
- Making checksums deterministic, this ensures that a compilation process will not be repeated unnecessarily.
- Added the dxWDL version number to the properties.
- Optimized queries that if a workflow/applet has already been compiled. This is done by limiting
the possible names of the data object we are searching for. We add a regular expression to the
query, bounding the legal names.

## 1.17
**Fixed**
- Setting recurse to false, when searching for existing applets and
  workflows. Fix contributed by Jeff Tratner.
- An optional output file, that does not exist, is returned as `None`.
For example, examine task `foo`.

```wdl
version 1.0
task foo {
    command {}
    output {
        File? f = "A.txt"
        Array[File?] fa = ["X.txt", "Y.txt"]
    }
}
```

Running it will return:
```
{
    "fa" : [null, null]
}
```
Note that `f` is missing. When passing through the DNAx
system, it is removed becaues it is optional and null.


## 1.16
**Fixed**
- [Bug 284](https://github.com/dnanexus/dxWDL/issues/284), dropping a default when
calling a subworkflow.

## 1.15

Improved find-data-objects queries, reducing the time to check if an
applet (or workflow) already exists on the platform. This is used when
deciding if an applet should be built, rebuilt, or archived.

To speed the query so it works on large projects with thousands of
applets and workflows, we limited the search to data objects generated
by dxWDL. These have a `dxWDL_checksum` property. This runs the risk
of missing cases where an applet name is already in use by a regular
dnanexus applet/workflow. We assume this is an unusual case. To fix this,
the existing applet can be moved or renamed.

## 1.14

**Fixed**
- Binary and decimal units are respected when specifying memory. For example, GB is 10<sup>9</sup> bytes, and GiB is 2<sup>30</sup> bytes. This follows the [memory spec](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#memory).

**Added**
- Initial support for the development WDL version, 1.1 (or 2.0). This does not include the directory type.

## 1.13
**Fixed**
- [Bug 274](https://github.com/dnanexus/dxWDL/issues/274), losing pipe symbols ('|') at the beginning of a line.

**Changed**
- Upgrade to Cromwell 44, with support for JSON-like values in meta sections

**Added**
- Ability to put a list of upstream projects into the [extras file](./doc/ExpertOptions.md#Setting-dnanexus-specific-attributes-for-tasks).

## 1.12
**Fixed**
- Tolerate platform applets/workflows with input/output specs that use non WDL types. For example, an array of applets.
- Bug when accessing call results that were not executed. For example, in `path_not_taken` the `compare` call is not made. This, incorrectly, causes an exception to be raised while evaluating `equality`.
```
version 1.0

workflow path_not_taken {
    if (false) {
        call compare
    }
    output {
        Boolean? equality = compare.equality
    }
}

task compare {
    command {}
    output {
        Boolean equality = true
    }
}
```

**Changed**
- Removed warning for variables that cannot be set from the inputs file. These are messages like this:
```
Argument unify.contig_shards, is not treated as an input, it cannot be set
Argument etl.shards, is not treated as an input, it cannot be set
Argument seq.iter_compare, is not treated as an input, it cannot be set
```
Top level call argument can now be set.

## 1.11
**Fixed**
- Default values specified for top-level calls using a JSON file.
- [Bug 272](https://github.com/dnanexus/dxWDL/issues/272)
- Runtime error when using an array of WDL structs

**Changed**
- Upgrade to Cromwell v43
- Replaced internal implementation of finding source lines for calls, with Cromwell/WOM implementation.

## 1.10
**Fixed**
- Handling of pair type in a JSON input file
- Allow overriding default values in calls from the input file. For example, in a workflow like:

```wdl
version 1.0

workflow override {
    call etl { input: a = 3 }
    output {
        Int result = etl.result
    }
}

task etl {
    input {
        Int a
        Int b = 10
    }
    command {}
    output {
        Int result = a + b
    }
}
```

We can now set b to a value other than 10, with an input file like this:

```
{
  "override.etl.b" : 5
}
```

**New**
- Support for compressed docker images (gzip)

**Changed**
- Upgrade to Cromwell v42


## 1.09
**Fixed**
- The `-p` flag was not respected
- Environment not passed correctly during compilation.
- Sorting structs by dependencies

## 1.08
**Fixed**
- [bug 259](https://github.com/dnanexus/dxWDL/issues/259). Unified the code for resolving platform paths.
- Error when downloading a file that resides in a container, rather than a project.
- Imports specified from the command line

**Changed**
- Merged `DxPath` and `DxBulkResolve` modules.

## 1.06
- Upgrade to Cromwell v41
- Mark auxiliary workflows and applets as hidden. This is a step in the direction of supporting copying of a workflow from one project to another.
- Added unit tests for the `WomValueAnalysis` module, which checks if WOM expressions are constant.
- Reducing reliance on the dxjava library, calling the DXAPI directly.
- Merged fork of the dxjava library back into the dx-toolkit.
- Use the dx-download-agent (dxda) instead of `dx download` in tasks.
- Fix for bug (https://github.com/dnanexus/dxWDL/issues/254)
- Fix for bug occurring when a `struct` is imported twice

## 1.05
- Got DxNI to work for DNAnexus apps, not just applets.
- Added links to applets that are referenced inside WDL fragments. This should allow, at some point, copying workflows between projects.
- Moved applet meta-information into the `details` field. This makes dx:applets more readable.
- Wrote unit-tests for the `WdlVarLinks` module.
- Batching file-describe, and file-resolve calls.

## 1.04
- Precalculate instance types. There are tasks that calculate the instance type they need in the runtime section. If the task call is made from a workflow fragment, we can calculate the instance type then and there. This allows launching directly in the correct instance type, instead of launching an additional job.
- Fixed bug in expression evaluation at runtime (https://github.com/dnanexus/dxWDL/issues/240)
- Improved unit-test coverage

## 1.03
- Bug fixes related to optional arguments (https://github.com/dnanexus/dxWDL/issues/235).

## 1.02
- Removed the instance-type database, and wom source code from the inputs of applets.
- Added the WDL source code to workflow and applet objects on the platform. It is stored in the details field, and
can be easily [retrieved](./doc/ExpertOptions.md#Getting-WDL-sources). It has been removed from the generated applet bash script.
- Initial support for the `struct` type
- Check that the reserved substring '___' is not used in the source WDL code. This sequence is used
  to translate dots ('.') into DNAx inputs and outputs. Dots are invalid symbols there.
- Bug fixes: https://github.com/dnanexus/dxWDL/issues/224, https://github.com/dnanexus/dxWDL/issues/227, https://github.com/dnanexus/dxWDL/issues/228

## 1.01
- Ensure that native docker uses the machine's hostname (i.e., the job ID) as
  the hostname, matching the previous behavior of dx-docker. This allows
  setting the job ID in error messages, helping debugging. Contributed by Jeff Tratner.

## 1.00
- Support for WDL version 1.0, as well as draft-2. There are two features
that are not yet supported: `struct`, and nested scatters. Work is ongoing
to address these omissions.

## 0.81.5
- Adding debug information to conversions from JSON to WDL variables. This helps
track down runtime problems with missing variables that a task is expecting. Previously,
we didn't know which variable was having a problem.

## 0.81.4
- Bug fix for complex cases where WDL files import each other.

## 0.81.3
- Bug fix for native docker. There was a problem when a docker image was using an internally defined
user; it didn't have the necessary permissions to create output files on the worker.

## 0.81.2
- Support NTLM proxies. If your organization is configured with an NTLM proxy,
you can use it like this:
```
$ export HTTP_PROXY_METHOD=ntlm
$ export HTTP_PROXY_DOMAIN = acme.com
$ export HTTP_PROXY = https://john_smith:welcome1@proxy.acme.com:8080
$ java -jar dxWDL.jar ...
```


## 0.81.1
- Supporting proxy configuration with user and password. For example:
```
$ export HTTPS_PROXY = https://john_smith:welcome1@proxy.acme.com:8080
```

## 0.81

- [Proxy configuration](./doc/ExpertOptions.md#Proxy-configurations). If your organization interposes a proxy between the internal machines and external hosts, you can set the environment variable `HTTP_PROXY` (or `HTTPS_PROXY`) to point to the proxy. The compiler will pass all of its dnanexus API calls through that proxy. For example, if you perform the following on the command line shell:

```bash
$ export HTTP_PROXY = proxy.acme.com:8080
$ java -jar dxWDL.jar ...
```

the compiler will route all requests through the machine `proxy.acme.com` on port `8080`.

## 0.80

- Native docker is now the default. If you still want to use [dx-docker](https://wiki.dnanexus.com/Developer-Tutorials/Using-Docker-Images), the `-useDxDocker` flag is available. In order to store a docker image on the platform, you can do `docker save`, and upload the tarball to a file. More details are provided in the [Expert options](./doc/ExpertOptions.md#Docker).

- The compiler emits a warning for partialy defined workflow outputs. For example, in workflow `foo`, the output `add.result` is partial, because it is not assigned to a variable. Partial definitions are discarded during compilation, hence the warning.

```wdl
workflow foo {
  call add { ... }
  output {
     add.result
  }
}
```

To avoid this problem, rewrite like this:

```wdl
workflow foo {
  call add { ... }
  output {
     Int r = add.result
  }
}
```

- Update test scripts for python3

## 0.79.1

- Version [docker-based runner script](https://github.com/dnanexus/dxWDL/tree/master/scripts/compiler_image/run-dxwdl-docker)
- Do not call sudo in runner script, in case system is set up not to require
  sudo to run Docker.
- Rename run script to `run-dxwdl-docker`

```
$ export DXWDL_VERSION=0.79.1
$ sudo run-dxwdl-docker compile /path/to/foo.wdl -project project-xxx
```

## 0.79
- Support [per task dx-attributes](./doc/ExpertOptions.md#Setting-dnanexus-specific-attributes-for-tasks).
- Report a warning for a non-empty runtime section in a native applet, instead of throwing an error. Note
that the WDL runtime section will be ignored, the native definitions will be used instead.
- Fix bug when using [spaces in output files](https://github.com/dnanexus/dxWDL/issues/181)
- Eliminate a job-describe API call from all tasks. This reduces overall platform load,
which is important in volume workflows.
- Support for [private docker registries](./doc/ExperOptions.md#Native-docker-and-private-registries)
- A [docker image](https://hub.docker.com/r/dnanexus/dxwdl) for
running the compiler without needing to install dependencies. You can
use the
[dxwdl](https://github.com/dnanexus/dxWDL/tree/master/scripts/compiler_image/dxwdl)
script to compile file `foo.wdl` like this:

```bash
$ runc.sh compile /path/to/foo.wdl -project project-xxxx
```

- Instructions for how to replace the built in
  [reorganize workflow outputs](./doc/ExpertOptions.md#Use-your-own-applet) applet,
  with your [own](./doc/ExpertOptions.md#Handling-intermediate-workflow-outputs).


## 0.78.1
- Support the `restartableEntryPoints` applet option in the `extras` file.

## 0.78
- Clone the dxWDL runtime asset to local project, to allow sub-jobs access to it.

## 0.77
- Improve user message when pretty printing an erroneous WDL file.
- New command line flag `--leaveWorkflowsOpen`, that leaves the toplevel workflow
open. This option is intended for power users, it allows modifying the workflow
after the compiler is done.
- Figure out the smallest instance from the *current* price list,
avoid the use of a hardcoded instance type. On the Azure cloud, an
applet failed because it tried using the hardcoded `mem1_ssd1_x4`
instance, which does not exist there (only on AWS).

## 0.76
- Handle using an asset that lives in another project by creating a local
record for it.

## 0.75
- Upgrade to Ubuntu 16.04
- Preparatory work for supporting WOM as a compilation intermediate representation.

## 0.74.1
- Removed inlining of tasks code into auxiliary applets. This was causing
a workflow to change when a task was modified, violating separate compilation,
and not allowing executable cloning to work.

## 0.74
- Improving test and release scripts.
- Adding printout for the dxWDL version to running applets and workflows
- Check and throw an exception if an asset is not in the current project. It
needs to be cloned.
- Supporting Amsterdam (azure:westeurope) and Berlin (aws:eu-central-1) regions.
- Fixed error in collecting results from an optional workflow branch.

## 0.72
- Put the project-wide reuse of applets under a special flag `projectWideReuse`.
- Improved the queries to find dx:executables on the target path.
- Improvements to the algorithm for splitting a WDL code block into parts.

## 0.71
- In an unlocked workflow, compile toplevel calls with no
subexpressions to dx stages. The
[expert options](./doc/ExpertOptions.md#toplevel-calls-compiled-as-stages)
page has a full description of the feature.
- Allow the compiler to reuse applets that have been archived. Such
applets are moved to a `.Archive` directory, and the creation date is
appended to the applet name, thereby modifying the applet name. The
name changes causes the search to fail. This was fixed by loosening the
search criterion.

## 0.70
- Upgrade to Cromwell 33.1
- Reuse applets and workflows inside a rpoject. The compiler now looks
for an applet/workflow with the correct name and checksum anywhere in
the project, not just in the target directory. This resolved
issue (https://github.com/dnanexus/dxWDL/issues/154).

## 0.69
- Support importing http URLs
- Simplified handling of imports and workflow decomposition
- Fixed issue (https://github.com/dnanexus/dxWDL/issues/148). This was a bug
occuring when three or more files with the same name were downloaded to
a task.
- Fixed issue (https://github.com/dnanexus/dxWDL/issues/146). This
occurred when (1) a workflow called a task, and (2) the task and the workflow had
an input with the same name but different types. For example, workflow
`w` calls task `PTAsays`, both use input `fruit`, but it has types
`Array[String]` and `String`.

```wdl
workflow w {
    Array[String] fruit = ["Banana", "Apple"]
    scatter (index in indices) {
        call PTAsays {
            input: fruit = fruit[index], y = " is good to eat"
        }
        call Add { input:  a = 2, b = 4 }
    }
}

task PTAsays {
    String fruit
    String y
     ...
}
```

## 0.68.1
- Fixing build and release scripts.
- Allowing non-priviliged users to find the dxWDL runtime asset in
the public repository.

## 0.68
- Throw an exception if a wrapper for a native platform call has
a non-empty runtime section.
- Use an SSD instance for the collect sub-jobs.
- Remove the runtime check for calling a task with missing values,
fix issue (https://github.com/dnanexus/dxWDL/issues/112). The check
is overly restrictive. The task could have a default, or, be able to
tolerate the missing argument.

## 0.67
- Color coding outputs, yellow for warnings, red for errors.
- Indenting information output in verbose mode
- Removing the use of `dx pwd`. The user needs to specify the
destination path on the command line. A way to avoid this, is to compile
from the dx-toolkit, with the upcoming `dx compile` command.


## 0.66.2
- Allow variable and task names to include the sub-string "last"
- Workflow fragment applets inherit access properties from the extras
file.

## 0.66.1
- Ignoring unknown runtime attributes in the extras file.

## 0.66
- Decomposing a workflow when there is a declaration after a call. For example,
workflow `foo` needs to be decomposed. The workflow fragment runner does not
handle dependencies, and cannot wait for the `add` call to complete.

```wdl
workflow foo {
    Array[Int] numbers

    scatter (i in numbers) {
        call add { input: a=i, b=1}
        Int m = add.result + 2
    }

    output {
        Array[Int] ms = m
    }
}
```

- Setting debug levels at runtime. The compiler flag `runtimeDebugLevel` can be set to 0, 1, or 2.
Level 2 is maximum verbosity, level 1 is the default, zero is minimal outputs.

- Upgrade to Cromwell version 32.

- Using headers to speed up the workflow decomposition step. The idea
is to represent a WDL file with header. When a file is imported, we
use the header, instead of pulling in the entire WDL code, including
its own imports.

- Support setting defaults in applets, not just workflows. This can be done with the `--defaults`
command line option, and a JSON file of WDL inputs.

## 0.65
- Optimization for the case of launching an instance where it is
calculated at runtime.

- Support dnanexus configuration options for tasks. Setting
the execution policy, timeout policies, and access
control can be achieved by specifying the default option in the
`default_taskdx_attributes` section of the `extras` file. For
example:

```json
{
  "default_task_dx_attributes" : {
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

- Improved error message for namespace validation. Details are no longer hidden when
the `-quiet` flag is set.

- Reduced logging verbosity at runtime. Disabled printing of directory structure when running
tasks, as the directories could be very large.

- Added support for calling native DNAx apps. The command

```
java -jar dxWDL.jar dxni -apps -o my_apps.wdl
```

instructs the compiler to search for all the apps you can call, and create WDL
tasks for them.


## 0.64
- Support for setting defaults for all task runtime attributes has been added.
This is similar to the Cromwell style. The `extras` command line flag takes a JSON file
as an argument. For example, if `taskAttrs.json` is this file:
```json
{
    "default_runtime_attributes" : {
      "docker" : "quay.io/encode-dcc/atac-seq-pipeline:v1"
    }
}
```

Then adding it to the compilation command line will add the `atac-seq` docker image to all
tasks by default.
```
java -jar dxWDL-0.44.jar compile test/files.wdl -defaults test/files_input.json -extras taskAttrs.json
```

- Reduced the number of auxiliary jobs launched when a task specifies the instance type dynamically.
A task can do this is by specifiying runtime attributes with expressions.
- Added the value iterated on in scatters

## 0.63
- Upgrade to cromwell-31 WDL/WOM library
- Report multiple validation errors in one step; do not throw an exception for
the first one and stop.
- Reuse more of Cromwell's stdlib implementation
- Close generated dx workflows


## 0.62.2
- Bug fix for case where optional optional types were generated. For example, `Int??`.

## 0.62

- Nested scatters and if blocks
- Support for missing arguments has been removed, the compiler will generate
an error in such cases.


## 0.61.1

- When a WDL workflow has an empty outputs section, no outputs will be generated.

## 0.61

- Decomposing a large block into a sub-workflow, and a call. For
example, workflow `foobar` has a complex scatter block, one that has
more than one call. It is broken down into the top-level workflow
(`foobar`), and a subworkflow (`foobar_add`) that encapsulates the
inner block.

```
workflow foobar {
  Array[Int] ax

  scatter (x in ax) {
    Int y = x + 4
    call add { input: a=y, b=x }
    Int base = add.result
    call mul { input: a=base, n=x}
  }
  output {
    Array[Int] result = mul.result
  }
}

task add {
  Int a
  Int b
  command {}
  output {
    Int result = a + b
  }
}

task mul {
  Int a
  Int n
  command {}
  output {
    Int result = a ** b
  }
}
```

Two pieces are together equivalent to the original workflow.
```
workflow foobar {
  Array[Int] ax
  scatter (x in ax) {
    Int y = x + 4
    call foobar_add { x=x, y=y }
  }
  output {
    Array[Int] result = foobar_add.mul_result
  }
}

workflow foobar_add {
  Int x
  Int y

  call add { input: a=y, b=x }
  Int base = add.result
  call mul { input: a=base, n=x}

  output {
     Int add_result = add.result
     Int out_base = base
     Int mul_result = mul.result
   }
 }
```

- A current limitation is that subworkflows, created in this way, give errors
for missing variables.
- The allowed syntax for workflow outputs has been tightened. An output
declaration must have a type and and a value. For example, this is legal:
```
output {
   Int add_result = add.result
}
```

but this is not:
```
output {
   add.result
}
```

## 0.60.2
- Minor bug fixes

## 0.60.1
- Bug fix release

## 0.60
- Split README into introduction, and advanced options. This should, hopefully, make
the top level text easier for the beginner.
- A workflow can call another workflow (subworkflow). Currently,
  the UI support for this feature is undergoing improvements.
- The search path for import can be extended by using the `--imports` command line argument.
This is useful when the source WDL files are spread across several directories.

## 0.59
- Improved handling of pricing list
- Upgraded to the new wdl4s library, that now lives in the cromwell
  repository under cromwell-wdl. The cromwell version is 30.2.

## 0.58
- Adhering to WDL spec: a declaration set to a constant is
treated as an input. In the example, `contamination` is compiled to
a workflow input with a default of `0.75`.

```
workflow w {
  Float contamination = 0.75
}
```

- Improving naming of workflow stages, and how their appear in the UI.
- Fixed regression in compilation of stand-alone applets.
- Fixed bug when pretty-printing applets that use <<<,>>> instead of {,}

## 0.57
- Fixed bug in setting workflow defaults from JSON file
- Improved behavior of the `sub` and `size` stdlib functions
- Improved naming convention for scatter/if workflow stages.

## 0.56
- Access to arguments in calls inside scatters. This
feature was lost while supporting locked-down workflows.

```
task Add {
    Int a
    Int b
    command {}
    output {
        Int result = a + b
    }
}

workflow w {
    scatter (i in [1, 10, 100]) {
        call Add { input: a=i }
    }
}
```

For example, in workflow `w`, argument `b` is missing from the `Add`
call. It is now possible to set it it from the command line with `dx
run w -iscatter_1.Add_b=3`. With an inputs file, the syntax is:
`w.Add.b = 3`.

- Fixed bad interaction with the UI, where default values were omitted.

## 0.55
- Improved support for conditionals and optional values
- Automated tests for both locked and regular workflows
- Various minor bug fixes

## 0.54
- More accurate detection of the IO classes in dx:applets.
- Improved handling of WDL constants. For example, multi word
constants such as `["1", "2", "4"]`, and `4 + 11` are recognized as
such. When possible, they are evaluated at runtime.
- Revert default compilation to a regular workflow. It is possible to
compile to a locked-down workflow by specifying `-locked` on the
command line. To specify command workflow inputs from the command line
use:
```
dx run foo -i0.N=19
```


## 0.53
- Not closing dx:workflow objects, this is too restrictive.


## 0.52
- Uplift the code to use DNAnexus workflow inputs and outputs.
Running a workflow now has slighly easier command line syntax.
For example, if workflow `foo` takes an integer
argument `N`, then running it through the CLI is done like
this:

```
dx run foo -iN=19
```

- Revamp the conversions between dx:file and wdl:file. This allows
specifying dx files in defaults.
- Initial support for non-empty WDL array types, for example `Array[String]+`.
- Improved handling of optional values

- An experimental new syntax for specifying inputs where the workflow
is underspecified. WDL allows leaving required call inputs unassigned, and
specifying them from the input file. For example, workflow `math`
calls task `add`, but does not specify argument `b`. It can then
be specified from the input file as follows: `{ "math.add.b" : 3}`.

```
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
       add.result
    }
}
```

The dx:workflow that is compiled from `math` can set this variable from
the command line as follows:
```
dx run math -iadd___b=5
```

This command line syntax may change in the future.

## 0.51
- Adding a `-quiet` flag to suppress warning and informational
  messages during compilation.
- Minor fixes

## 0.50
- Tasks that have an empty command section do not download files eagerly. For example,
if a task only looks at `size(file)`, it will not download the file at all.

```
task fileSize {
    File in_file

    command {}
    output {
        Float num_bytes = size(in_file)
    }
}
```

- Support docker images that reside on the platform.


## 0.49
- Removed limitation on maximal hourly price per instance. Some bioinformatics pipelines
regularly require heavy weight instances.
- Fixed several WDL typing issues
- Improving the internal representation of JSON objects given as input.

## 0.48
- Support string interpolation in workflow expressions
- Handle scatters where the called tasks return file arrays, and non-native
platform types. This is done by spawning a subjob to collect all the outputs
and bundle them.
- Azure us-west region is supported

## 0.47
- Support calling native DNAx applets from a WDL workflow. A helper utility
  is `dxni` (*Dx Native Interface*), it creates task wrappers for existing
  dx:applets.

## 0.46
- Allow applet/workflow inputs that are optional, and have a default.
- More friendly command line interface

## 0.45
- Default workflow inputs. The `--defaults` command line argument
embeds key-value pairs as workflow defaults. They can be overridden
at runtime if necessary.

## 0.44
- Use hashes instead of files for non-native dx types
- Do not use the help field in applet input/output arguments to carry
  WDL typing information.
- Fold small text files into the applet bash script. This speeds up the 'dx build'
  command, and avoids uploading them as separate platform files.
- Replace 'dx build' with direct API calls, additional compilation speedup.

## 0.43
- Speed up the creation of a dx:workflow. With one API call create
all the workflow stages, and add set the checksum property.
- Speeding up applet construction

## 0.42
- WDL objects: experimental
- Type coercion
- Faster compilation when applets already exist. Using bulk describe API, instead
  of looking up each applet separately.
- Adding checksum to dx:workflow objects, recompile only if it has not changed
- Specified behavior of input/output files in doc/Internals.md

## 0.41
- Minor fixes to streaming

## 0.40
- Upgrade to wdl4s version 0.15.
- Enable several compiler warnings, removing unused imports
- Bug fixes:
  * Wrong dx project name when using multiple shells, and
    different projects
  * Better handling of errors from background 'dx cat' processes.

## 0.39
- Streaming works with tasks that use docker

## 0.38
- Minor fixes, and better testing, for file streaming

## 0.37
- File download streaming

## 0.36
- Conditionals

## 0.35
- Support white space in destination argument
- Allow providing the dx instance type directly in the runtime attributes.

## 0.34
- Topological sorting of workflow. By default, check for circular
dependencies, and abort in case there are cycles. Optionally, sort the
calls to avoid forward references. This is useful for WDL scripts
that were not written by hand.
- Create an `output_section` applet for a compiled workflow. If
`--reorg` is specified on the command line, the applet also
reorganizes the output folder; it moves all intermediate results to
subfolder `intermediate`. Reorganization requires `CONTRIBUTE` level access,
which the user needs to have.
- Multi region support.

## 0.33
- Initial support for WDL pairs and maps
- Upgrade to wdl4s version 0.13
- Improve parsing for memory specifications
- Optionals can be passed as inputs and outputs from tasks.

## 0.32
- Support archive and force flags a la `dx build`. Applets and workflows
  are *not* deleted by default, the --force flag must be provided.
- If there are insufficient permissions to get the instance price list, we
  have a reasonable fallback option.
- Added namespace concept to intermediate representation.
- WDL pretty printer now retains output section.

## 0.31
- The compiler is packed into a single jar file
- Fixed glob bug

## 0.30
Bug fix release.

Corrected the checksum algorithm, to recurse through the directory
tree.

## 0.29
Improved support for scatters
  * Support an expression in a scatter collection.
  * Compile declarations before a scatter into the scatter
  applet. This optimization folds two applets into one.
  * Use wdl4s native name resolution to calculate scatter
  collection type.
  * Added documentation to doc/IR.md describing how scatters
  are compiled.

## 0.28
- Support for ragged file arrays
- Correctly handle an empty workflow output section
- Upgrade to wdl4s version 0.12
- Improvements to regression tests

## 0.27
- Support the import directive. This allows spreading definitions
  across multiple files.
- Using native wdl4s data structures, instead of pretty printing, for
  internal representation.

## 0.26
- Support passing unbound variables to calls inside scatters
- Implemented glob and size function for files
- Correctly handling empty arrays in task input/output
- Improved scatter linking. Got rid of runtime calls to get dx:applet descriptions.

## 0.25
- Allow a docker image to be specified as a parameter
- Use Travis containers for continuous integration
- Adding tests for instance type choice at runtime

## 0.24
- Support resource requirements (memory/disk) calculated from runtime variables.

## 0.23
- Rebuild an applet, only if the source code has changed. For example,
if an applet's WDL code changed, it will be rebuilt in the next compilation.

## 0.22
- Files are lazily downloaded in auxiliary applets, such as scatters. In
  tasks/applets, they are always downloaded.

## 0.21
- Adding linking information to applets that call other applets. This
ensures the correct applet-ids are invoked. No lookup by name
is performed at runtime.
- Minor bug fixes

## 0.20
- Separate compilation of workflows and tasks. A task is compiled to
single dx:applet, and a workflow is compiled to auxiliary applets and
a dx:workflow.
- WDL files containing only tasks, and no workflow, are supported
- Many improvements all around

## 0.13
- Upgrade to wdl4s v0.10
- Improving version number handling
- Minor improvements in debugging printout and error reporting

## 0.12
- The GATK pipeline compiles and runs (result not validated yet)
- Bug fixes:
  * Correct order of scatter output arrays
  * Improved handling of optionals

## 0.11
- Better compilation error messages
- Version checking
- Compiler verbose mode
- Renamed initial workflow stage to Common
