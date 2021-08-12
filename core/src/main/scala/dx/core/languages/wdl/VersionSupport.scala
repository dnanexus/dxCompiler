package dx.core.languages.wdl

import java.nio.file.Path
import dx.api.DxApi
import dx.core.languages.wdl.WdlUtils.parseSource
import dx.util.{Bindings, FileNode, FileSourceResolver, Logger, StringFileNode}
import spray.json._
import wdlTools.generators.code.WdlGenerator
import wdlTools.syntax.{Parsers, WdlParser, WdlVersion}
import wdlTools.types.{Section, TypeCheckingRegime, WdlTypes, TypedAbstractSyntax => TAT}

case class WdlOptions(regime: TypeCheckingRegime.TypeCheckingRegime) {
  def toJson: JsValue = {
    JsObject(
        "regime" -> JsString(regime.toString)
    )
  }
}

object WdlOptions {
  val default: WdlOptions = WdlOptions(regime = TypeCheckingRegime.Moderate)

  def fromJson(jsValue: JsValue): WdlOptions = {
    jsValue match {
      case JsNull => WdlOptions.default
      case JsObject(fields) =>
        fields.get("regime") match {
          case Some(JsString(value)) => WdlOptions(TypeCheckingRegime.withNameIgnoreCase(value))
          case None                  => WdlOptions.default
          case other =>
            throw new Exception(s"invalue 'regime' value ${other}")
        }
      case other =>
        throw new Exception(s"invalid WdlOptions ${other}")
    }
  }
}

case class VersionSupport(version: WdlVersion,
                          wdlOptions: WdlOptions = WdlOptions.default,
                          fileResolver: FileSourceResolver = FileSourceResolver.get,
                          dxApi: DxApi = DxApi.get,
                          logger: Logger = Logger.get,
                          wdlParser: Option[WdlParser] = None) {
  private lazy val parser = wdlParser.getOrElse(
      Parsers(followImports = true, fileResolver = fileResolver, logger = logger).getParser(version)
  )

  lazy val generatedVersion: WdlVersion = {
    version match {
      case WdlVersion.Draft_2 =>
        logger.warning("Upgrading draft-2 input to version 1.0")
        WdlVersion.V1
      case _ => version
    }
  }

  lazy val codeGenerator: WdlGenerator = WdlGenerator(Some(generatedVersion))

  def parse(
      sourceCode: FileNode,
      wdlOptions: WdlOptions = WdlOptions.default
  ): (TAT.Document, Bindings[String, WdlTypes.T_Struct]) = {
    WdlUtils.parseAndCheckSource(sourceCode, parser, wdlOptions, fileResolver, logger)
  }

  def parse(src: String): (TAT.Document, Bindings[String, WdlTypes.T_Struct]) = {
    parse(StringFileNode(src))
  }

  def parse(path: Path): (TAT.Document, Bindings[String, WdlTypes.T_Struct]) = {
    parse(fileResolver.fromPath(path))
  }

  def parseExpression(exprStr: String,
                      bindings: Bindings[String, WdlTypes.T],
                      docSource: FileNode,
                      section: Section.Section = Section.Other): TAT.Expr = {
    WdlUtils.parseExpr(exprStr, version, parser, docSource, bindings, section)
  }

  def validateWdlCode(wdlWfSource: String): Unit = {
    val (tDoc, _) = parse(wdlWfSource)
    // Check that this is the correct language version
    if (tDoc.version.value != version) {
      throw new Exception(
          s"document has wrong version ${tDoc.version.value}, should be ${version}"
      )
    }
  }

  def generateDocument(doc: TAT.Document): String = {
    val sourceString = codeGenerator.generateDocument(doc).mkString("\n")
    Logger.get.ignore(
        WdlUtils.parseAndCheckSourceString(sourceString,
                                           doc.source.toString,
                                           wdlOptions,
                                           fileResolver,
                                           logger)
    )
    sourceString
  }

  def generateElement(element: TAT.Element): String = {
    val sourceString = codeGenerator.generateElement(element).mkString("\n")
    // add the version statement so we can try to parse it
    val standAloneString = s"version ${generatedVersion.name}\n\n${sourceString}"
    // we only do parsing, not type checking here, since the element may
    // not be stand-alone
    val parser = Parsers.default.getParser(generatedVersion)
    Logger.get.ignore(parseSource(StringFileNode(standAloneString), parser))
    standAloneString
  }

}

object VersionSupport {
  def fromSource(
      sourceCode: FileNode,
      wdlOptions: WdlOptions = WdlOptions.default,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      dxApi: DxApi = DxApi.get,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, WdlTypes.T_Struct], VersionSupport) = {
    val parser = Parsers(followImports = true, fileResolver = fileResolver, logger = logger)
      .getParser(sourceCode)
    val (doc, typeAliases) =
      WdlUtils.parseAndCheckSource(sourceCode, parser, wdlOptions, fileResolver, logger)
    val versionSupport =
      VersionSupport(doc.version.value, wdlOptions, fileResolver, dxApi, logger, Some(parser))
    (doc, typeAliases, versionSupport)
  }

  def fromSourceFile(
      sourceFile: Path,
      wdlOptions: WdlOptions = WdlOptions.default,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      dxApi: DxApi = DxApi.get,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, WdlTypes.T_Struct], VersionSupport) = {
    fromSource(fileResolver.fromPath(sourceFile), wdlOptions, fileResolver, dxApi, logger)
  }

  def fromSourceString(
      sourceCode: String,
      wdlOptions: WdlOptions = WdlOptions.default,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      dxApi: DxApi = DxApi.get,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, WdlTypes.T_Struct], VersionSupport) = {
    fromSource(StringFileNode(sourceCode), wdlOptions, fileResolver, dxApi, logger)
  }
}
