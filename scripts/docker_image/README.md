# dxCompiler Docker Image

## Build

- Specify the dxCompiler version to use [here](./docker.env).
- Run [docker_build](./docker_build.sh)

## Run

- Login to DNAnexus with `dx login` and select a project.
- Specify the dxCompiler version to use [here](./docker.env).
- Run [docker_run](./docker_run.sh) with dxCompiler args as args.

```bash
./docker_run.sh version
```

```bash
./docker_run.sh -help  
```
