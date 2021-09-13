# Developing dxCompiler

## Setting up your development environment

* Install JDK 8
    * On mac with [homebrew](https://brew.sh/) installed:
      ```
      $ brew tap AdoptOpenJDK/openjdk
      $ brew cask install adoptopenjdk8
      # Use java_home to find the location of JAVA_HOME to set
      $ /usr/libexec/java_home -V
      # We recommend adding the following line to ~/.zshrc
      $ export JAVA_HOME=/Library/Java/...
      ```
    * On Linux (assuming Ubuntu 16.04)
      ```
      $ sudo apt install openjdk-8-jre-headless
      ```
    * Scala will compile with JDK 11, but the JDK on DNAnexus worker instances is JDK 8 and will not be able to run a JAR file with classes compiled by a later version of Java
* Install [sbt](https://www.scala-sbt.org/), which also installs Scala. Sbt is a make-like utility that works with the ```scala``` language.
    * On MacOS: `brew install sbt`
    * On Linux:
      ```
      $ wget www.scala-lang.org/files/archive/scala-2.12.1.deb
      $ sudo dpkg -i scala-2.12.1.deb
      $ echo "deb https://dl.bintray.com/sbt/debian /" | sudo tee -a /etc/apt/sources.list.d/sbt.list
      $ sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
      $ sudo apt-get update
      $ sudo apt-get install sbt
      ```
    * Running sbt for the first time takes several minutes, because it downloads all required packages.
* We also recommend to install [Metals](https://scalameta.org/metals/), which enables better integration with your IDE
    * For VSCode, install the "Scala (Metals)" and "Scala Syntax (official)" plugins
* You will also need Python installed (we recommend version 3.6+) and dx-toolkit (`pip install dxpy`)
* You will also need to create a GitHub personal access token (this is required by the sbt-github-packages plugin). In GitHub settings, go to "Developer settings > Personal access token" and create a new token with "write:packages" and "read:packages" scopes only. Then, export the `GITHUB_TOKEN` environment variable with this token as the value. For example, in your `.profile`:
  ```bash
  export GITHUB_TOKEN=<your personal access token>
  ```
    * On macOS, you may also want to add this token into your global environment so it is visible to your IDE:
  ```bash
  launchctl setenv GITHUB_TOKEN $GITHUB_TOKEN
  ```

## Getting the source code

* Clone or fork the [dxCompiler repository](https://github.com/dnanexus/dxCompiler) (depending on whether you have commit permissions)
* Checkout an existing dxCompiler branch or create a new branch (e.g. feat/42-my-feature)
* Add pre-commit hooks:
    * Create/edit a file .git/hooks/pre-commit
    * Add the following lines
      ```bash
      #!/bin/bash
      check=$(sbt scalafmtCheckAll)
      if [[ "$?" -ne "0" ]]; then
        echo "Reformatting; please commit again"
        sbt scalafmtAll
        exit $check
      fi
      ```
    * Make the file executable (e.g. `chmod +x .git/hooks/pre-commit`)
    * This hook runs the code formatter before committing code. You can run this command manually, but it is easiest just to have it run automatically.

## Developing in a Docker container

A Dockerfile with all the dependencies to build and test dxCompiler is available [here](docker/Dockerfile). To build an image from it and run a Docker container, run from the [docker](docker/) directory:

```
make
```

See below on how to run unit and integration tests. To recompile dxCompiler with your updates, run an integration test, as described below.

> Always make sure you push your changes to the remote github repo before destroying your Docker container.

## Adding new code

1. Checkout the `develop` branch.
2. Create a new branch with your changes. Name it something meaningful, like `APPS-123-download-bug`.
3. If the current snapshot version matches the release version, increment the snapshot version.
- For example, if the current release is `1.0.0` and the current snapshot version is `1.0.0-SNAPSHOT`, increment the snapshot version to `1.0.1-SNAPSHOT`.
- You can use a script to update the version simultaneously in all of the sub-packages: `scripts/update_version.md <version>`
4. Make your changes. Test locally using `sbt test`.
5. Update the release notes under the top-most header (which should be "in develop").
6. If the current snapshot version only differs from the release version by a patch, and you added any new functionality (vs just fixing a bug), increment the minor version instead.
- For example, when you first created the branch you set the version to `1.0.1-SNAPSHOT`, but then you realized you needed to add a new function to the public API, change the version to `1.1.0-SNAPSHOT`.
7. When you are done, create a pull request against the `develop` branch.

While developing, make sure you do the following:

* Follow the style guidelines (below).
* Always write unit tests for any new code you add, and update tests for any code you modify.
    * Unit tests should assert, and not print to the console
    * WDL test files belong in the top directory `test`

### Style guidelines

* We use [scalafmt style](https://scalameta.org/scalafmt/) with a few modifications. You don't need to worry so much about code style since you will use the automatic formatter on your code before it is committed.
* Readability is more important than efficiency or concision
    - write more/slower code if it makes the code more readable
    - avoid using abbreviations in variable names unless you know they will be obvious to everyone
    - add comments to clarify any non-obvious code
* Avoid using more complex features, e.g. reflection

## Using sbt

sbt is the build system for Scala (similar to a Makefile). The following are the main commands you will use while developing.

### Compiling the code

Scala (like Java) is a compiled language. To compile, run:

```
$ sbt compile
```

If there are errors in your code, the compiler will fail with (hopefully useful) error messages.

### Cleaning up artifacts & building

Generate a staging token via the web UI and login with `dx login --staging --token <token>`.

Run [scripts/clean_build.sh](scripts/clean_build.sh) to clean up existing artifacts (locally and on staging) and build new dxCompiler artifacts.

### Running unit tests

You should always run the unit tests after every successful compile. Generally, you want to run `sbt testQuick`, which only runs the tests that failed previously, as well as the tests for any code you've modified since the last time you ran the tests. However, the first time you checkout the code (to make sure your development environment is set up correctly) and then right before you push any changes to the repository, you should run the full test suite using `sbt test`.

You need to have a DNAnexus account and be logged into DNAnexus via the command line before you can run the tests (`dx login`).

### Running the integration tests

Integration tests actually build and run apps/workflows on DNAnexus. These tests take much longer to run than the unit tests, and so typically you only run them before submitting a pull request. You can also submit a PR and trigger the integration tests via Github Actions.

### Running integration tests on GitHub

You can run run integration tests after submitting a PR. By default the integration tests pipeline is skipped and only runs when the `integration` label is addded to the PR and in subsequent commit pushes. If you want to push more changes and temporarily skip these tests, remove the label.

The results will be available in the [Actions](https://github.com/dnanexus/dxCompiler/actions) tab. Ideally set the label only before requesting a review so that we don't incur too high costs from running the jobs at each push.

Note that only DNAnexus developers can set up a label on a PR so let us know when you'd like to request a review and we'll start them for you.

### Running integration tests locally

First, you need to an account on the DNAnexus staging environment, and you need to be added to the projects that have been setup for testing. Next, log into the DNAnexus staging environment using dx-toolkit: `dx login --staging`. Note that very often your login will timeout while the integration tests are running unless you are actively using the platform in another session, and this will cause the tests to fail. To avoid, this, generate a token via the web UI and use that token to log in on the command line: `dx login --staging --token <token>`.

Follow the "Cleaning up artifacts & building" instructions above.

Finally, run the integration tests. From the root dxCompiler directory, run `./scripts/run_tests.py`. You can run a specific test/suite by adding the `--test <name>` option.

Note that the test script does a lot of things for you. If for some reason you want to run them manually, here is what happens:

The dxCompiler and dxExecutor* JAR files are built and staged in the root dxCompiler directory. To do this manually, run `sbt assembly`, then move the JAR files from the `applet_resources` folder to the root dxCompiler folder, e.g. `mv applet_resources/dxCompiler.jar ./dxCompiler-X.Y.Z.jar`.

### Running a subset of tests locally

You can also execute a subset of tests, for example to run medium set of tests:

```
./scripts/run_tests.py --test M
```

It will also generate a runtime asset in the test project and a dxCompiler*.jar on your local computer so you can use them for further manual testing of the compilation and execution.

It's also possible to specify one test to run from the [/test](/test) directory and invoke it by name, for example:

```
./scripts/run_tests.py --test add3
```

Check the test runner script `--help` for more options.

### Test data on platform

Any files that tests rely on should be stored in `dxCompiler_playground:/test_data/` on `staging`.

### Other sbt tips

#### Cache

sbt keeps the cache of downloaded jar files in `${HOME}/.ivy2/cache`. For example, the WDL jar files are under `${HOME}/.ivy2/cache/org.broadinstitute`. In case of problems with cached jars, you can remove this directory recursively. This will make WDL download all dependencies (again).

## Releasing a new version

### Releasing using Github Actions

dxCompiler can be released from Github. The release pipeline (optionally) runs large integration tests, builds the release on staging, runs multi-region tests on staging (one test per region), builds on production, and creates a Docker image, which is pushed to DockerHub.

1. Checkout the develop branch (either HEAD or the specific commit you want to release)
2. Create a release branch named with the version number, e.g. `release-2.4.2`
3. Update the version numbers in application.conf files by removing `-SNAPSHOT`
    - Run `scripts/update_version.md <version>`
    - If you want to update the versions manually, there are 5 of them:
        * [compiler](https://github.com/dnanexus/dxCompiler/blob/main/compiler/src/main/resources/application.conf)
        * [core](https://github.com/dnanexus/dxCompiler/blob/main/core/src/main/resources/application.conf)
        * [executorCommon](https://github.com/dnanexus/dxCompiler/blob/main/executorCommon/src/main/resources/application.conf)
        * [executorWdl](https://github.com/dnanexus/dxCompiler/blob/main/executorWdl/src/main/resources/application.conf)
        * [executorCwl](https://github.com/dnanexus/dxCompiler/blob/main/executorCwl/src/main/resources/application.conf)
4. Update the [Release Notes](https://github.com/dnanexus/dxCompiler/blob/main/RELEASE_NOTES.md)
    - Change the top header from "in develop" to "<version> (<date>)"
5. Update versions of libraries as needed in [build.sbt](https://github.com/dnanexus/dxCompiler/blob/main/build.sbt).
6. Push the release branch to GitHub.
7. Run the release pipeline:
    1. Go to `Actions` > `dxCompiler Release (Staging and Prod)` and click `Run workflow` on the right side.
    2. Make sure the `release-xxx` branch is selected (default setting).
    3. Once finished, the pipeline will create a draft release page on GitHub.
8. Test in customer projects: build and run the applet [test-dxcompiler-in-project](/test-dxcompiler-in-project) in customer projects that are shared with `org-dnanexus_apps_customer_testers`.
9. Publish the draft [release](https://github.com/dnanexus/dxCompiler/releases). The compressed source code (in `zip` and `tar.gz`) will be added to the release page automatically.

If you encounter any additional issues while creating the release, you will need to make the fixes in `develop` and then merge them into the release branch.

To complete the release, open a PR to merge the release branch into main. You can then delete the release branch.

Following the release, it is also nice to bump the SNAPSHOT versions in the `develop` branch.

### Releasing manually

This should only be done if you want to create a debug release for internal testing (and even then, you can follow the automated process above and just not publish the draft release).

1. Follow steps 1-4 above
2. Make sure the Unit and Integration tests and [customer tests](/test-dxcompiler-in-project) pass.
3. Clean your `dx` environment because you'll be using limited-power tokens to run the release script. Do not mix them with your regular user token.
    ```
    dx clearenv
    ```
4. Build new externally visible release
    ```
    ./scripts/build_all_releases.sh --staging-token XXX --production-token YYY --docker-user UUU --docker-password WWW
    ```
   this will take a while. It builds the release on staging, runs multi-region tests on staging (one test per region), builds on production, and creates an easy to use Docker image, which is pushed to DockerHub.

## Enabling dxCompiler in new DNAnexus regions

If dxCompiler needs to be enabled in a new DNAnexus region the following should be done in the staging & production environments:

* dxCompiler requires WDL and CWL assets to be stored in a public "dxCompiler\_region" project. Please use these [scripts](/scripts/dxCompiler_projects_utils) to create a new project and add a "PUBLIC" invitee to the project.
* The release script that tests and releases dxCompiler in all regions needs to be [updated](/scripts/build_release.py#L33) (one line change).
* The app used for copying the assets to different regions during a release ([app-dxwdl_copy](/scripts/dxcompiler_copy)) needs to be enabled in the new region (please update `regionalOptions`, `whatsNew`, and increment the `version` of the app).
