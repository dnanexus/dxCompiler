# dxCompiler Docker Image

Steps to build and run a Docker image with dxCompiler, as an alternative to installing Java and running the dxCompiler.jar locally.

## Build Docker

1. Specify the dxCompiler version to use [here](./docker.env).
2. Run [docker_build](./docker_build.sh)

## Run Docker

1. Login to DNAnexus with `dx login` and select a project.
2. Run [docker_run](./docker_run.sh) with dxCompiler args as args.

```bash
./docker_run.sh -help  
```

- For compiling a workflow, the local working directory is mounted to the container.

```bash
# General example
<path to docker_run.sh> compile <path to WDL source>

# From dxCompiler repo directory
./scripts/docker_image/docker_run.sh compile test/single_tasks/add3.wdl
```
