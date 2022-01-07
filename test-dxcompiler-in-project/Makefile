# TODO: find a way to download from draft release URLs; the below works for final releases.
# For draft releases, save jar manually to resources/home/dnanexus/dxCompiler.jar
VERSION=2.8.3
DXCOMPILER_PRERELEASE_PATH=project-G7GFgjQ0X386z5ZF24vvyFzp:dxCompiler-${VERSION}.jar
DXCOMPILER_RELEASE_URL=https://github.com/dnanexus/dxCompiler/releases/download/latest/dxCompiler-${VERSION}.jar

all: download-prerelease

no-op:
	echo "Using dxCompiler.jar in resource directory"

download-prerelease:
	mkdir -p resources/home/dnanexus
	dx download ${DXCOMPILER_PRERELEASE_PATH} -o resources/home/dnanexus/dxCompiler.jar

download-release:
	mkdir -p resources/home/dnanexus
	wget ${DXCOMPILER_JAR_URL} -O resources/home/dnanexus/dxCompiler.jar
