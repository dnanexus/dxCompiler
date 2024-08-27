# dxCompiler Docker Image

## Build Docker

- Specify the dxCompiler version to use [here](./docker.env).
- Run [docker_build](./docker_build.sh)

## Run Docker

- Login to DNAnexus with `dx login` and select a project.
- Specify the dxCompiler version to use [here](./docker.env).
- Run [docker_run](./docker_run.sh) with dxCompiler args as args.

```bash
./docker_run.sh -help  
```

- For compiling a workflow, the local working directory is mounted to the container

```bash
# General example
<path to docker_run.sh> compile <path to WDL source>

# From dxCompiler repo directory
./scripts/docker_image/docker_run.sh compile test/single_tasks/add3.wdl
```
