# An image that allows compiling a WDL program. It saves you the
# requirement to install java.

FROM ubuntu:16.04
ARG VERSION

# Install java-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk-headless

# dxCompiler
COPY dxCompiler-$VERSION.jar /dxCompiler.jar

ENTRYPOINT ["java", "-jar", "/dxCompiler.jar"]
