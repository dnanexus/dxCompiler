FROM amazoncorretto:11
ARG VERSION
ARG BASE_DIR

COPY ${BASE_DIR}/dxCompiler-$VERSION.jar /dxCompiler.jar

ENTRYPOINT ["java", "-jar", "/dxCompiler.jar"]
