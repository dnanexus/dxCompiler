# dxCint

## Motivation
`dxCint` (**dxC**ompiler **int**egration) is a new framework for integration tests for dxCompiler. It allows to easily 
add new tests, run tests and extend the framework with new test classes/types. It also includes a set of unit tests for 
self testing.

## Installing
1. Clone this repo
2. `cd dxcint`
3. `poetry install` (recommended, because will be installed from `poetry.lock` with pinned dependency versions). Or 
`pip install .`

## CLI methods
#### Adding tests
See [Adding tests](#adding-tests-1) section
```bash
Usage: dxcint add [OPTIONS] DXC_REPOSITORY_ROOT SUITE CATEGORY

  Add tests to the suite
  Positional Arguments:
      DXC_REPOSITORY_ROOT: A root directory of a dxCompiler repository. Should contain build.sbt.
      SUITE: Test suite name. Usually a team-defined group of tests to be run in each CI/CD step. Check README for the
              CLI command to extract available suites.
      CATEGORY: Test category name. Usually a team-defined category that reflects a type of the test.
              Check dxcint/config/{SUITE_FILE}.json, where the keys are category names.

Options:
  -d, --directory TEXT  Directory containing workflow fixtures to be added.
                        DEFAULT: current WD
  -x, --extension TEXT  Extension of the test workflow script files. Allowed
                        are 'cwl', 'cwl.json', 'wdl'. DEFAULT: wdl
```

#### Running tests
See [Running tests](#running-tests-1) section
```bash
Usage: dxcint integration [OPTIONS] DXC_REPOSITORY_ROOT

  Run integration test or test suite
  Positional Arguments:
      DXC_REPOSITORY_ROOT: A root directory of a dxCompiler repository. Should contain build.sbt.

Options:
  -t, --test_name TEXT  Name of a single test or suite of tests. If not
                        provided - builds just the core compiler library (e.g.
                        dxCompiler-VERSION.jar)

```

#### Available suites
Command to print a hash of available suites
```bash
Usage: dxcint suites [OPTIONS]
```

Returns the structure:
```bash
Available suites are:
{"Suite name": "config file name in the dxcint/config where test names are registered within categories"}
```

## Running tests
To run the tests just call 
```bash
dxcint integration ${ROOT_DIR_OF_DXCOMPILER_REPO} --test_name ${TEST_NAME}
```
where `${TEST_NAME}` is the name of a single test ot a test suite. Available test suites can be checked by `dxcint suites`, 
see [Available suites](#available-suites) for more info.


## Adding tests
1. Create a workflow/task/tool script file intended as a test fixture. **Convention**: file name and the workflow name must be identical. 
Workflow names and task/tool names must be unique and not repeating any existing names throughout the existing tests. 
**Naming hint**: if working on a JIRA ticket, include it in the workflow and task names (e.g. `workflow apps_123_buggy_bug`)
2. Create - if appropriate - input and result fixtures using test name as a prefix and `_input.json` or `_results.json` 
suffixes respectively. Example - `apps_123_buggy_bug_input.json`. Place these files in the same directory as the test 
workflow file (see #1).
3. Run `dxcint add` with appropriate arguments. This will add your test to the test suite config file (one of the files 
in `dxcint/config`) according to specified test suite and category. Respective file fixtures (workflow, inputs and 
results) will be copied to `dxcint/resources`.

Tests can be added manually to the suite by placing them in the appropriate directory in `dxcint/resources` and adding 
the names in appropriate config file. **But this is not recommended.**

## Test categories
All tests are `-locked` by default, except `unlocked_expected_output`

### Basic categories
|   test category    | implementation class name  |                                  description                                   |
|:------------------:|:--------------------------:|:------------------------------------------------------------------------------:|
|  registered_test   |       RegisteredTest       |                       always fails, parent of all tests                        |
| analysis_finished  |      AnalysisFinished      |                 `_validate` checks that the analysis succeeded                 |
|  expected_failure  |      ExpectedFailure       |                  `_validate` checks that the analysis failed                   |
|  expected_output   |       ExpectedOutput       |  `_validate` checks that the expected result fixture matches analysis outputs  |


### Complex categories
Complex categories are obtained by adding [mixins](#mixins) to the [basic categories](#basic-categories)

|             test category              |     implementation class name      |                                                            comments                                                             |
|:--------------------------------------:|:----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------:|
|        expected_failure_message        |       ExpectedFailureMessage       | ExpectedFailure + ResultsMixin, `_validate` checks that the analysis failed with an error message specified in results['error'] |
|        extras_analysis_finished        |       ExtrasAnalysisFinished       |                                                 AnalysisFinished + ExtrasMixin                                                  |
|         extras_expected_output         |        ExtrasExpectedOutput        |                                                  ExpectedOutput + ExtrasMixin                                                   |
|        unlocked_expected_output        |       UnlockedExpectedOutput       |                                   ExpectedOutput + UnlockedMixin, _extract_outputs overriden                                    |
|       manifest_analysis_finished       |      ManifestAnalysisFinished      |                                                AnalysisFinished + ManifestMixin                                                 |
|         reorg_expected_output          |        ReorgExpectedOutput         |                                                   ExpectedOutput + ReorgMixin                                                   |
| static_pinned_instance_expected_output | StaticPinnedInstanceExpectedOutput |                                     ExpectedOutput + StaticOnlyMixin + PinnedInstanceMixin                                      |
|         extern_expected_output         |        ExternExpectedOutput        |                  before compiling, creates an applet interface from additional applets in dxcint/dependencies                   |
|       app_extern_expected_output       |      AppExternExpectedOutput       |                                      before compiling, creates a task interface from apps                                       |


#### Mixins

Mixins are designed so that many can be used to extend a category class.

|        Name         |                                                     purpose                                                     |
|:-------------------:|:---------------------------------------------------------------------------------------------------------------:|
| PinnedInstanceMixin |                                        `instance_type` is set to `None`                                         |
|     ExtrasMixin     |                            appends to compile flags `--extras testname_extras.json`                             |
|    UnlockedMixin    |                                      removes `-locked` from compiler flags                                      |
|    ManifestMixin    |                                    appends to compile flags `-useManifests`                                     |
|     ReorgMixin      |                                        appends to compile flags `-reorg`                                        |
|  ResultsTestMixin   |              loads result fixture from `testname_results.json` and stores it in `results` property              |
|   StaticOnlyMixin   | appends to compile flags `-instanceTypeSelection static` preventing random selection of dynamic/static instance |


## Extending the framework

### Adding new test types
New test types can be easily added to the framework. First, one needs to add a new test category (e.g. 'expected_failure'). 
Then, if this category needs a new behavior when compiling, running and/or validating the test, a new subclass of 
`RegisteredTest` needs to be implemented. There might be up to 3 methods that might need implementation.

To do that, use [**mixins**](#mixins) (classes designed for multiple inheritance) to add common functionality to your 
new test type, put them first in the derived class definition. E.g. `class NewTestClass(Mixin1, Mixin2, AnalysisFinished)`. 

1. **Compilation**. If a test workflow in this new test category should be compiled with a specific combination of dxCompiler flags, 
the property `exec_id` should be overridden/implemented with passing those compiler arguments to the 
`_compile_executable()` method. 
2. **Execution**. To change executable run behavior - the methods `_run_executable()` and `_run_executable_inner()` can be overridden.
3. **Validation**. For test validation, each subclass has to have an implementation of `_validate()` method.


## Test suite architecture
The following architecture is implemented to support the integration tests for dxCompiler.
1. **CLI** orchestrates the flow 
2. **Context**: immutable class maintaining the passed arguments and values which dxcint was initiated with. 
3. **TestDiscovery** scans for config files, registers the tests by creating an array of RegisteredTest objects. 
4. **Dependency** provides the necessary dependencies for building dxCompiler assets on the platform.
5. **Terraform** takes the list of **RegisteredTest**s and prepares the DNAnexus platform for running them.
6. **RegisteredTest** has methods to compile test, run test and validate the results.
7. **Messenger** class to communicate between RegisteredTest and the respective analysis/job on the platform.


## Adding a dependency 
Create a new dependency config file to add a new dependency for your tests. Dependencies are packaged in the assets 
specific for each DSL (domain specific language - WDL, CWL, etc.) supported by dxCompiler. Dependency config is a `.json` file with basic information about a 
dependency package. Developers have to create one config file for each dependency they need to be included into the 
language specific assets for dxCompiler.  
**!!! Please consult the IT and compliance teams before adding new dependencies**
#### Structure
**Required fields**:
1. "name": str. Name of a dependency. E.g. "awscli"
2. "languages": Array[str]. Languages that this dependency is built for. E.g. ["wdl", "cwl"]
3. "version": str. Version of the software to be used.

Required fields that define type of the dependency (implemented as separate subclasses of `Dependency`)
* "package_manager": str. Also, all necessary tags should be added as [here](https://documentation.dnanexus.com/developer/apps/execution-environment#external-utilities).
* "source_link": str. Link where to get the binary executable, if it was not downloaded before.  
The "package_manager" attribute has a precedence over "source_link". If "package_manager" is not provided, "source_link" 
will be used to download the binary executable of the dependency. If neither is used dxCint will assume the dependency 
is installed on the workers by default and nothing will happen.


