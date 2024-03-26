# TODO: find a way to download from draft release URLs; the below works for final releases.
# For draft releases, save jar manually to resources/home/dnanexus/dxCompiler.jar
DXCOMPILER_JAR_URL=https://github.com/dnanexus/dxCompiler/releases/download/2.11.5/dxCompiler-2.11.5.jar

all: no-op

no-op:
	echo "Using dxCompiler.jar in resource directory"

download-latest:
	mkdir -p resources/home/dnanexus
	wget ${DXCOMPILER_JAR_URL} -O resources/home/dnanexus/dxCompiler.jar
