package dx.core.ir

import dx.AppInternalException
import dx.api.{DxApi, DxExecution, DxFile, DxFileDescCache, DxUtils, DxWorkflowStage}
import dx.core.Constants
import dx.core.ir.Value._
import spray.json._
import dx.util.{Enum, FileSourceResolver, Logger}

object IORef extends Enum {
  type IORef = Value
  val Input, Output = Value
}

/**
  * A union of all the different ways of building a value from JSON passed
  * by the platform. A complex value is a WDL values that does not map to
  * a native dx:type. Such values may also have files embedded in them.
  * For example:
  * - Ragged file array:  Array\[Array\[File\]\]
  * - Object with file elements
  * - Map of files:     Map[String, File]
  * A complex value is implemented as a json structure, and an array of
  * all the files it references.
  */
sealed trait ParameterLink {
  val dxType: Type

  /**
    * Copy this ParameterLink, replacing dxType with its optional equivalent.
    * @return
    */
  def makeOptional: ParameterLink
}
case class ParameterLinkValue(jsv: JsValue, dxType: Type) extends ParameterLink {
  def makeOptional: ParameterLinkValue = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

case class ParameterLinkStage(dxStage: DxWorkflowStage,
                              ioRef: IORef.Value,
                              varName: String,
                              dxType: Type)
    extends ParameterLink {
  def makeOptional: ParameterLinkStage = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

case class ParameterLinkWorkflowInput(varName: String, dxType: Type) extends ParameterLink {
  def makeOptional: ParameterLinkWorkflowInput = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

case class ParameterLinkExec(dxExecution: DxExecution, varName: String, dxType: Type)
    extends ParameterLink {
  def makeOptional: ParameterLinkExec = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

object ParameterLink {
  // Key used to wrap a complex value in JSON.

  val FlatFilesSuffix = "___dxfiles"
  val WorkflowInputFieldKey = "workflowInputField"
}

case class ParameterLinkSerializer(fileResolver: FileSourceResolver = FileSourceResolver.get,
                                   dxApi: DxApi = DxApi.get,
                                   pathsAsObjects: Boolean = false) {

  /**
    * Serialize a complex value into a JSON value. The value could potentially point
    * to many files. The assumption is that files are already in DNAnexus format,
    * so not requiring upload/download or any special conversion.
    * @param t the type
    * @param v the value
    * @return
    */
  private def serialize(t: Type, v: Value): JsValue = {
    if (Type.isNestedOptional(t)) {
      Logger.error(s"""|jsFromWdlValue
                       |    type=${t}
                       |    val=${v}
                       |""".stripMargin)
      throw new Exception("a nested optional type/value")
    }

    def handler(irValue: Value, irType: Type): Either[Value, JsValue] = {
      (irType, irValue) match {
        case (_, VString(s)) if s.length > Constants.StringLengthLimit =>
          throw new AppInternalException(
              s"string is longer than ${Constants.StringLengthLimit}"
          )
        case (Type.TMulti.Any, f: VFile) =>
          Right(
              ValueSerde.wrapValue(ValueSerde.serializePath(f, Some(fileResolver), pathsAsObjects))
          )
        case (_, f: VFile) =>
          Right(ValueSerde.serializePath(f, Some(fileResolver), pathsAsObjects))
        case (Type.TFile, VString(path)) =>
          Right(ValueSerde.serializePath(VFile(path), Some(fileResolver)))
        case (Type.TMulti.Any, d: DirectoryValue) =>
          Right(
              ValueSerde.wrapValue(ValueSerde.serializePath(d, Some(fileResolver), pathsAsObjects))
          )
        case (_, d: DirectoryValue) =>
          Right(ValueSerde.serializePath(d, Some(fileResolver), pathsAsObjects))
        case (Type.TDirectory, VString(path)) if Value.isDxFolderUri(path) =>
          Right(ValueSerde.serializePath(VFolder(path), Some(fileResolver), pathsAsObjects))
        case _ => Left(irValue)
      }
    }
    ValueSerde.serializeWithType(v, t, Some(handler))
  }

  /**
    * Create a link from a WDL value.
    *
    * @param t the WDL type
    * @param v the WDL value
    * @return
    */
  def createLink(t: Type, v: Value): ParameterLink = {
    ParameterLinkValue(serialize(t, v), t)
  }

  def createConstantField(value: Value,
                          bindName: String,
                          encodeDots: Boolean = true): (String, JsValue) = {
    val bindEncName =
      if (encodeDots) {
        Parameter.encodeName(bindName)
      } else {
        bindName
      }
    (bindEncName, ValueSerde.serialize(value))
  }

  def serializeSimpleLink(link: ParameterLink): JsValue = {
    link match {
      case ParameterLinkValue(jsLinkvalue, _) => jsLinkvalue
      case ParameterLinkStage(dxStage, ioRef, varEncName, _) =>
        ioRef match {
          case IORef.Input =>
            dxStage.getInputReference(Parameter.encodeName(varEncName))
          case IORef.Output =>
            dxStage.getOutputReference(Parameter.encodeName(varEncName))
        }
      case ParameterLinkWorkflowInput(varEncName, _) =>
        JsObject(
            DxUtils.DxLinkKey -> JsObject(
                ParameterLink.WorkflowInputFieldKey -> JsString(
                    Parameter.encodeName(varEncName)
                )
            )
        )
      case ParameterLinkExec(dxJob, varEncName, _) =>
        DxUtils.dxExecutionToEbor(dxJob, Parameter.encodeName(varEncName))
    }
  }

  // create input/output fields that bind the variable name [bindName] to
  // this WdlVar
  def createFieldsFromLink(link: ParameterLink,
                           bindName: String,
                           encodeDots: Boolean = true): Vector[(String, JsValue)] = {
    val encodedName =
      if (encodeDots) {
        Parameter.encodeName(bindName)
      } else {
        bindName
      }
    if (Type.isNative(link.dxType, !pathsAsObjects)) {
      // Types that are supported natively in DX
      Vector((encodedName, serializeSimpleLink(link)))
    } else {
      // Complex type requiring two fields: a JSON structure, and a flat array of files.
      val fileArrayName = s"${encodedName}${ParameterLink.FlatFilesSuffix}"
      val mapValue = link match {
        case ParameterLinkValue(jsLinkValue, _) =>
          // files that are embedded in the structure
          val jsFiles = JsArray(DxFile.findFiles(dxApi, jsLinkValue).map(_.asJson))
          // Dx allows hashes as an input/output type. If the JSON value is
          // not a hash (JsObject), we need to add an outer layer to it.
          val jsLink = JsObject(Parameter.ComplexValueKey -> jsLinkValue)
          Map(encodedName -> jsLink, fileArrayName -> jsFiles)
        case ParameterLinkStage(dxStage, ioRef, varName, _) =>
          val varFileArrayName = s"${varName}${ParameterLink.FlatFilesSuffix}"
          ioRef match {
            case IORef.Input =>
              Map(
                  encodedName -> dxStage.getInputReference(varName),
                  fileArrayName -> dxStage.getInputReference(varFileArrayName)
              )
            case IORef.Output =>
              Map(
                  encodedName -> dxStage.getOutputReference(varName),
                  fileArrayName -> dxStage.getOutputReference(varFileArrayName)
              )
          }
        case ParameterLinkWorkflowInput(varName, _) =>
          val varFileArrayName = s"${varName}${ParameterLink.FlatFilesSuffix}"
          Map(
              encodedName ->
                JsObject(
                    DxUtils.DxLinkKey -> JsObject(
                        ParameterLink.WorkflowInputFieldKey -> JsString(varName)
                    )
                ),
              fileArrayName ->
                JsObject(
                    DxUtils.DxLinkKey -> JsObject(
                        ParameterLink.WorkflowInputFieldKey -> JsString(varFileArrayName)
                    )
                )
          )
        case ParameterLinkExec(dxJob, varName, _) =>
          val varFileArrayName = s"${varName}${ParameterLink.FlatFilesSuffix}"
          Map(
              encodedName -> DxUtils
                .dxExecutionToEbor(dxJob, Parameter.encodeName(varName)),
              fileArrayName -> DxUtils
                .dxExecutionToEbor(dxJob, Parameter.encodeName(varFileArrayName))
          )
      }
      mapValue.toVector
    }
  }

  def createFields(bindName: String,
                   t: Type,
                   v: Value,
                   encodeDots: Boolean = true): Vector[(String, JsValue)] = {
    createFieldsFromLink(createLink(t, v), bindName, encodeDots)
  }

  def createFieldsFromMap(values: Map[String, (Type, Value)],
                          encodeDots: Boolean = true): Map[String, JsValue] = {
    values.flatMap {
      case (name, (t, v)) => createFields(name, t, v, encodeDots)
    }
  }
}

case class ParameterLinkDeserializer(dxFileDescCache: DxFileDescCache, dxApi: DxApi = DxApi.get) {
  private def unwrapComplex(jsv: JsValue): JsValue = {
    jsv match {
      case JsObject(fields) if fields.contains(Parameter.ComplexValueKey) =>
        // unpack the hash with which complex JSON values are wrapped in dnanexus.
        fields(Parameter.ComplexValueKey)
      case _ => jsv
    }
  }

  def deserializeInput(jsv: JsValue): Value = {
    ValueSerde.deserialize(unwrapComplex(jsv),
                           dxApi = dxApi,
                           dxFileDescCache = Some(dxFileDescCache))
  }

  def deserializeInputMap(inputs: Map[String, JsValue]): Map[String, Value] = {
    inputs.map {
      case (name, jsv) => name -> deserializeInput(jsv)
    }
  }

  def deserializeInputWithType(
      jsv: JsValue,
      t: Type,
      name: String,
      handler: Option[(JsValue, Type) => Either[JsValue, Value]] = None
  ): Value = {
    ValueSerde.deserializeWithType(unwrapComplex(jsv),
                                   t,
                                   name,
                                   handler,
                                   dxApi,
                                   Some(dxFileDescCache))
  }
}
