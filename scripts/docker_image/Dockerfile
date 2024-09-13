FROM amazoncorretto:11
ARG VERSION

COPY dxCompiler-$VERSION.jar /dxCompiler.jar

ENTRYPOINT ["java", "-jar", "/dxCompiler.jar"]
