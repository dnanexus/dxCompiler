FROM amazoncorretto:11
ARG VERSION
ARG BASE_DIR

COPY ${BASE_DIR}/dxCompiler-$VERSION.jar /dxCompiler.jar

# Copy the entrypoint script
COPY docker_image_entry_test.sh /docker_image_entry_test.sh
RUN chmod +x /docker_image_entry_test.sh

ENTRYPOINT ["/docker_image_entry_test.sh"]
