# How to set up dxCompiler in DNAnexus mini dev instance

Every applet compiled from WDL/CWL using dxCompiler has a reference to a runtime executor asset, which is required to properly execute the applet.
The asset is stored in a public `dxCompiler` project (one per region) in the `/releases/\<version\>` folder. Each version of dxCompiler contains a config
file that stores a json with a mapping of region projects and folder paths where the asset of the corresponding version is stored (the dxCompiler doesn't
store asset record IDs).

In order to move dxCompiler to a new environment (other than staging or prod), for example mini environment, the asset for the dxCompiler version you want to work with needs to be downloaded
from the DNAnexus platform, uploaded to your environment, and stored in a location where the dxCompiler executor expects it, as shown below.

## Upload the runtime executor asset to your environment

The following example assumes you will be working in the aws:us-east-1 region and uses version 2.10.0 as the version.

1. Download the asset from the DNAnexus platform (either prod or staging):
```
# Set environment variables
DXCOMPILER_VERSION=2.10.0
LANGUAGE=WDL
PROJECT_NAME=dxCompiler

RECORD=${PROJECT_NAME}:/releases/$DXCOMPILER_VERSION/dx${LANGUAGE}rt
ASSET_FILE=$(dx describe ${RECORD} | grep "Outgoing links" | awk '{print $3}')
dx download $ASSET_FILE
```
 
2. Create a project with the same name, a specific folder in it, and upload the asset to it:

```
# Set the same environment variables
DXCOMPILER_VERSION=2.10.0
LANGUAGE=WDL
PROJECT_NAME=dxCompiler

# Create project and folder
DXCOMPILER_PROJECT_ID=$(dx new project --region aws:us-east-1 ${PROJECT_NAME} --brief)
dx mkdir -p ${DXCOMPILER_PROJECT_ID}:/releases/$DXCOMPILER_VERSION

# Upload the asset
ASSET_FILE_ID=$(dx upload - --brief --visibility hidden --destination ${PROJECT_NAME}:/releases/$DXCOMPILER_VERSION/dx${LANGUAGE}rt.tar.gz)
 
# Create a record pointing to the asset
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
```     

## Download the compiler to your local machine

In order to compile workflows in your mini env, download the dxCompiler with the same version you downloaded & uploaded the asset, from:

```
DXCOMPILER_VERSION=2.10.0
wget https://github.com/dnanexus/dxCompiler/releases/download/$DXCOMPILER_VERSION/dxCompiler-$DXCOMPILER_VERSION.jar
```
