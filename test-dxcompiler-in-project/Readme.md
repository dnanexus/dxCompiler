# test-dxcompiler-in-project Readme

Applet `test-dxcompiler-in-project` can be built and deployed in customer-provided projects to facilitate testing on their WDL workflows.

## Test Project Setup

- The [example_project directory](./example_project) shows the recommended structure for a test project.
- `test/test_wdl.sh` should be a script that downloads the WDL files onto the worker, compiles the workflows, and runs them.
- It is recommended to have 1 `workflow_#` folder per test workflow containing the WDL and auxiliary files necessary to compile the workflow.
- It is recommended to direct compiled artifacts and workflow outputs to `test_out/workflow_#`.
- Input data files should be kept elsewhere; they do not need to be downloaded for compilation.

## Test Project Access

- Customer projects should be shared with `org-dnanexus_apps_customer_testers`.

## Building & Running the Applet

- Upload regional dxCompiler assets to the platform and stage dxCompiler.jar in the root `dxCompiler` folder. This can be done bu running the initial stage of a release.
- `dx select` the project for testing.
- `dx build` and `dx run` this applet.
