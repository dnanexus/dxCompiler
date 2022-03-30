# How to set up dxCompiler in DNAnexus mini dev instance

Every applet compiled from WDL/CWL using dxCompiler has a reference to a runtime executor asset, which is required to properly execute the applet.

The asset is stored in a public `dxCompiler` project (one per region) in the `/releases/<version>` folder. Each version of dxCompiler contains a config file that stores a json with a mapping of region projects and folder paths where the asset of the corresponding version is stored (the dxCompiler doesn't store asset record IDs).

In order to move dxCompiler to a new environment (other than staging or prod), for example mini environment, the asset for the dxCompiler version you want to work with needs to be downloaded from the DNAnexus platform, uploaded to your environment, and stored in a location where the dxCompiler executor expects it, as shown below.

## In DNAnexus staging or prod environment - download the runtime executor asset

The following example assumes you will be working in the aws:us-east-1 region and uses version 2.10.0 as the version.

1. First specify the version of dxCompiler you want to work with:

```
export DXCOMPILER_VERSION=2.10.0
```

2. Download the asset from the DNAnexus platform (either prod or staging):
```
#!/usr/bin/env bash

set -e

# Set environment variables
LANGUAGE=WDL
PROJECT_NAME=dxCompiler

RECORD=${PROJECT_NAME}:/releases/$DXCOMPILER_VERSION/dx${LANGUAGE}rt
ASSET_FILE=$(dx describe ${RECORD} | grep "Outgoing links" | awk '{print $3}')
echo "Downloading dxCompiler $DXCOMPILER_VERSION asset"
dx download $ASSET_FILE
echo "Done downloading $ASSET_FILE."
```
 
## In your DNAnexus mini env - upload the runtime executor asset to a specific location

1. Set the version
```
export DXCOMPILER_VERSION=2.10.0
```

2. Create a project with the same name as on prod/staging, a specific folder in it, and upload the asset to it

```
#!/usr/bin/env bash

set -e

# Set the same environment
LANGUAGE=WDL
PROJECT_NAME=dxCompiler

# Create project and folder
echo "Creating the project for dxCompiler $DXCOMPILER_VERSION asset"
DXCOMPILER_PROJECT_ID=$(dx new project --region aws:us-east-1 ${PROJECT_NAME} --brief)
dx mkdir -p ${DXCOMPILER_PROJECT_ID}:/releases/$DXCOMPILER_VERSION

# Upload the asset
echo "Uplading the asset dx${LANGUAGE}rt.tar.gz..""
ASSET_FILE_ID=$(dx upload dx${LANGUAGE}rt.tar.gz --brief --visibility hidden --destination ${PROJECT_NAME}:/releases/$DXCOMPILER_VERSION/dx${LANGUAGE}rt.tar.gz)

# Create a record pointing to the asset
echo "Creating the record.."
dx new record \
     --type AssetBundle \
     --details '{"archiveFileId": {"$dnanexus_link": "'"${ASSET_FILE_ID}"'"}}' \
     --property title="dx${LANGUAGE} asset" \
     --property description="Prerequisites for running ${LANGUAGE} workflows compiled to the platform" \
     --property version="${DXCOMPILER_VERSION}" \
     --property distribution="Ubuntu" \
     --property release="20.04" \
     --property runSpecVersion="0" \
     --close \
     ${PROJECT_NAME}:/releases/${DXCOMPILER_VERSION}/dx${LANGUAGE}rt

echo "Done."
```     

## Download the compiler to your local machine

In order to compile workflows in your mini env, download the dxCompiler **with the same version** you downloaded & uploaded the asset, from:

```
DXCOMPILER_VERSION=2.10.0
wget https://github.com/dnanexus/dxCompiler/releases/download/$DXCOMPILER_VERSION/dxCompiler-$DXCOMPILER_VERSION.jar
```

You can test your dxCompiler setup by compiling a small test WDL, for example the [add3 task](https://github.com/dnanexus/dxCompiler/blob/develop/test/single_tasks/add3.wdl) or the [echo_pairs workflow](https://github.com/dnanexus/dxCompiler/blob/develop/test/input_file/echo_pairs.wdl). For a quick intro on how to compile workflows with dxCompiler, please read the [Readme](https://github.com/dnanexus/dxCompiler).

