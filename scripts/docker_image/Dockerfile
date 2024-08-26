# dxCompiler in a Docker image. It saves you the requirement of installing Java.

FROM amazoncorretto:11
ARG VERSION

# dxCompiler
COPY dxCompiler-$VERSION.jar /dxCompiler.jar

ENTRYPOINT ["java", "-jar", "/dxCompiler.jar"]
