# An image that allows compiling a WDL program. It saves you the
# requirement to install java.

FROM ubuntu:16.04
ARG VERSION

# Install java-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk-headless

# dxCompiler
COPY dxCompiler-$VERSION.jar /dxCompiler.jar

# Copy the entrypoint script
COPY docker_image_entry.sh /docker_image_entry.sh
RUN chmod +x /docker_image_entry.sh

ENTRYPOINT ["/docker_image_entry.sh"]
