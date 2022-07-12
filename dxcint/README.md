# dxCint

## Motivation
`dxCint` (**dxC**ompiler **int**egration) is a new framework for integration tests for dxCompiler. It allows to easily 
add new tests, run tests and extend the framework with new test classes/types. It also includes a set of unit tests for 
self testing.

## Installing
1. Clone this repo
2. `cd dxcint`
3. `potry install` (recommended, because will be installed from `poetry.lock` with pinned dependency versions). Or 
`pip install .`

## CLI methods
#### Adding tests
See [Adding tests](#adding-tests) section
```commandline
Usage: dxcint add [OPTIONS] DXC_REPOSITORY_ROOT

Options:
  -d, --directory TEXT  Extension of the test workflow script files. Allowed
                        are 'cwl', 'cwl.json', 'wdl'
  -x, --extension TEXT  Extension of the test workflow script files. Allowed
                        are 'cwl', 'cwl.json', 'wdl'
  -s, --suite TEXT      Test suite name. Usually a team-defined group of tests
                        to be run in each CI/CD step
  -c, --category TEXT   Test category name. Usually a team-defined category
                        that reflects a type of the test
```

#### Running tests
See [Running tests](#running-tests) section
```commandline
Usage: dxcint integration [OPTIONS] DXC_REPOSITORY_ROOT

Options:
  -t, --test_name TEXT  Name of a single test or suite of tests. If not
                        provided - builds just the core compiler library (e.g.
                        dxCompiler-VERSION.jar)
```

#### Available suites
Command to print a hash of available suites with the structure:
```python
{"Suite name": "config file name in the dxcint/config where test names are registered"}
```
```commandline
Usage: dxcint suites [OPTIONS]

Options:
  --help  Show this message and exit.
```

## Running tests
To run the tests just call `dxcint integration`. Call with `--test_name` argument and provide a test name or a test 
suite to run a specific test or a suite of tests. Available test suites can be checked by `dxcint suites`, see 
[Available suites](#available-suites) for more info.


## Adding tests
1. Create a script file with the workflow intended as a test. **Convention**: file name and workflow name must be identical. 
Workflow name and task names must be unique and not repeating any existing names among the existing tests. 
**Naming hint**: if working on a JIRA ticket, include it in the workflow and task names (e.g. `workflow apps_123_buggy_bug`)
2. Create - if appropriate - input and result fixtures using test name as a prefix and `_input.json` or `_results.json` 
suffixes respectively. Example - `apps_123_buggy_bug_input.json`. Place these files in the same directory as the test 
workflow file (see #1).
3. Run `dxcint add` with appropriate arguments. This will add your test to the test suite config file (one of the files 
in `dxcint/config`) according to specified test suite and category. Respective file fixtures (workflow, inputs and 
results) will be copied to `dxcint/resources`.

Tests can be added manually to the suite by placing them in the appropriate directory in `dxcint/resources` and adding 
the names in appropriate config file.

## Extending the framework
### Adding new test types
New test types can be easily added to the framework. First, one needs to add a new test category (e.g. 'expected failure'). 
Then, if this category needs a new behavior when compiling, running and/or validating the test, a new subclass of 
`RegisteredTest` needs to be implemented. There might be up to 3 methods that might need implementation.

1. **Compilation**. If a test workflows in this new test category should be compiled with a specific combination of dxCompiler flags, 
the property `exec_id` should be overridden/implemented with passing those compiler arguments to the 
`_compile_executable()` method. 
2. **Execution**. To change executable run behavior - the methods `_run_executable()` and `_run_executable_inner()` can be overridden.
3. **Validation**. For test validation, each subclass has to have an implementation of `_validate()` method.


## Test suite architecture
The following architecture will be implemented to support integration tests for the compiler.
1. **CLI** orchestrates the flow 
2. **TestDiscovery** scans for config files, registers the tests by creating an array of RegisteredTest objects. 
3. **Dependency** class sets up necessary environment with immutable state for **RegisteredTest**’s.
4. **Terraform** takes the list of **RegisteredTest**’s and prepares the platform for running.
5. **RegisteredTest** has methods to compile test, run test and validate results.


## Dependency config 
Dependency config is a `.json` file with basic information about the dependency package. Developers have to create one 
for each dependency they need to be included into language specific assets for dxCompiler.  
**!!!Please consult the IT and compliance teams before adding new dependencies**
#### Structure
**Required fields**:
1. "name": str. Name of a dependency. E.g. “awscli”
2. "languages": Array[str]. Languages that this dependency is built for. E.g. [“wdl”, “cwl”]
3. "version": str. Version of the software to be used.

Required fields that define type of the dependency (implemented as separate subclasses of `Dependency`)
* "package_manager": str. Also, all necessary tags should be added as [here](https://documentation.dnanexus.com/developer/apps/execution-environment#external-utilities)
* "source_link": str. Link where to get the binary executable, if it was not downloaded before
The "package_manager" attribute has a precedence over "source_link". If "package_manager" is not provided, "source_link" 
will be used to download the binary executable of the dependency. If neither is used dxCint will assume the dependency 
is installed on the workers by default and nothing will happen.


