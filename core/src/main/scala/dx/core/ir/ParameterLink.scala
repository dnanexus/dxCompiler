package dx.core.ir

import dx.AppInternalException
import dx.api.{DxApi, DxExecution, DxFile, DxFileDescCache, DxUtils, DxWorkflowStage}
import dx.core.Constants
import dx.core.ir.Value._
import dx.util.protocols.DxFolderSource
import dx.util.{Enum, FileSourceResolver}
import spray.json._

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

object ParameterLink {
  // key used to wrap a complex value in JSON
  val WorkflowInputFieldKey = "workflowInputField"

  def prettyFormat(links: Map[DxName, ParameterLink], indent: Int = 2): String = {
    links
      .map { case (dxName, link) => s"${dxName} -> ${link}" }
      .mkString(s"\n${" " * indent}")
  }
}

case class ParameterLinkValue(jsv: JsValue, dxType: Type) extends ParameterLink {
  def makeOptional: ParameterLinkValue = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

/**
  * A ParameterLink that is a reference to another parameter.
  */
sealed trait ParameterLinkReference extends ParameterLink {

  /**
    * The type of the parameter source, which is usually, but not always,
    * the same as `dxType`. For example, a workflow input with type `String`
    * may be connected to a stage input with type `Any`.
    */
  val sourceType: Type
}

case class ParameterLinkStage(dxStage: DxWorkflowStage,
                              ioRef: IORef.Value,
                              dxName: DxName,
                              dxType: Type,
                              sourceType: Type)
    extends ParameterLinkReference {
  def makeOptional: ParameterLinkStage = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

case class ParameterLinkWorkflowInput(dxName: DxName, dxType: Type, sourceType: Type)
    extends ParameterLinkReference {
  def makeOptional: ParameterLinkWorkflowInput = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

case class ParameterLinkExec(dxExecution: DxExecution,
                             dxName: DxName,
                             dxType: Type,
                             sourceType: Type)
    extends ParameterLinkReference {
  def makeOptional: ParameterLinkExec = {
    copy(dxType = Type.ensureOptional(dxType))
  }
}

object ParameterLinkExec {
  def create(dxExecution: DxExecution, dxName: DxName, dxType: Type): ParameterLinkExec = {
    ParameterLinkExec(dxExecution, dxName, dxType, dxType)
  }
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
      throw new Exception(s"Trying to serialize a nested optional type ${t} for value ${v}")
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
        case (Type.TDirectory, VString(path)) if DxFolderSource.isDxFolderUri(path) =>
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

  def serializeSimpleLink(link: ParameterLink): JsValue = {
    link match {
      case ParameterLinkValue(jsLinkvalue, _) => jsLinkvalue
      case ParameterLinkStage(dxStage, ioRef, dxName, _, _) =>
        ioRef match {
          case IORef.Input =>
            dxStage.getInputReference(dxName.encoded)
          case IORef.Output =>
            dxStage.getOutputReference(dxName.encoded)
        }
      case ParameterLinkWorkflowInput(dxName, _, _) =>
        JsObject(
            DxUtils.DxLinkKey -> JsObject(
                ParameterLink.WorkflowInputFieldKey -> JsString(dxName.encoded)
            )
        )
      case ParameterLinkExec(dxJob, dxName, _, _) =>
        DxUtils.dxExecutionToEbor(dxJob, dxName.encoded)
    }
  }

  // create input/output fields that bind the variable name [bindName] to this parameter
  def createFieldsFromLink(link: ParameterLink, bindName: DxName): Vector[(DxName, JsValue)] = {
    val destTypeIsNative = Type.isNative(link.dxType, !pathsAsObjects)
    link match {
      case ref: ParameterLinkReference
          if destTypeIsNative && !Type.isNative(ref.sourceType, !pathsAsObjects) =>
        throw new Exception(s"cannot link non-native type to native type: ${ref}")
      case _ if destTypeIsNative =>
        // links with types that are supported natively in DX
        Vector((bindName, serializeSimpleLink(link)))
      case ref: ParameterLinkReference if Type.isNative(ref.sourceType, !pathsAsObjects) =>
        // References where the source type is native but the dest type is not.
        // For example, linking a `String` workflow input to a CWL `Any` stage input -
        // the link is still "simple" because there is no `___dxfiles` parameter associated
        // with the workflow input.
        Vector((bindName, serializeSimpleLink(link)))
      case _ =>
        // Complex type requiring two fields: a JSON structure, and a flat array of files.
        val fileArrayName = bindName.addSuffix(Constants.FlatFilesSuffix)
        val mapValue = link match {
          case ParameterLinkValue(jsLinkValue, _) =>
            // files that are embedded in the structure
            val jsFiles = JsArray(DxFile.findFiles(dxApi, jsLinkValue).map(_.asJson))
            // Dx allows hashes as an input/output type. If the JSON value is
            // not a hash (JsObject), we need to add an outer layer to it.
            val jsLink = JsObject(Constants.ComplexValueKey -> jsLinkValue)
            Map(bindName -> jsLink, fileArrayName -> jsFiles)
          case ParameterLinkStage(dxStage, ioRef, dxName, _, _) =>
            val varFileArrayName = dxName.addSuffix(Constants.FlatFilesSuffix)
            ioRef match {
              case IORef.Input =>
                Map(
                    bindName -> dxStage.getInputReference(dxName.encoded),
                    fileArrayName -> dxStage.getInputReference(varFileArrayName.encoded)
                )
              case IORef.Output =>
                Map(
                    bindName -> dxStage.getOutputReference(dxName.encoded),
                    fileArrayName -> dxStage.getOutputReference(varFileArrayName.encoded)
                )
            }
          case ParameterLinkWorkflowInput(dxName, _, _) =>
            val varFileArrayName = dxName.addSuffix(Constants.FlatFilesSuffix)
            Map(
                bindName ->
                  JsObject(
                      DxUtils.DxLinkKey -> JsObject(
                          ParameterLink.WorkflowInputFieldKey -> JsString(dxName.encoded)
                      )
                  ),
                fileArrayName ->
                  JsObject(
                      DxUtils.DxLinkKey -> JsObject(
                          ParameterLink.WorkflowInputFieldKey -> JsString(varFileArrayName.encoded)
                      )
                  )
            )
          case ParameterLinkExec(dxJob, dxName, _, _) =>
            val varFileArrayName = dxName.addSuffix(Constants.FlatFilesSuffix)
            Map(
                bindName -> DxUtils.dxExecutionToEbor(dxJob, dxName.encoded),
                fileArrayName -> DxUtils.dxExecutionToEbor(dxJob, varFileArrayName.encoded)
            )
        }
        mapValue.toVector
    }
  }

  def createFields(bindName: DxName,
                   t: Type,
                   v: Value
                   //encodeName: Boolean = true
  ): Vector[(DxName, JsValue)] = {
    createFieldsFromLink(createLink(t, v), bindName)
  }

  def createFieldsFromMap(
      values: Map[DxName, (Type, Value)]
      //encodeName: Boolean = true
  ): Map[DxName, JsValue] = {
    values.flatMap {
      case (name, (t, v)) => createFields(name, t, v)
    }
  }
}

case class ParameterLinkDeserializer(dxFileDescCache: DxFileDescCache, dxApi: DxApi = DxApi.get) {
  private def unwrapComplex(jsv: JsValue): JsValue = {
    jsv match {
      case JsObject(fields) if fields.contains(Constants.ComplexValueKey) =>
        // unpack the hash with which complex JSON values are wrapped in dnanexus.
        fields(Constants.ComplexValueKey)
      case _ => jsv
    }
  }

  def deserializeInput(jsv: JsValue): Value = {
    ValueSerde.deserialize(unwrapComplex(jsv),
                           dxApi = dxApi,
                           dxFileDescCache = Some(dxFileDescCache))
  }

  def deserializeInputMap(inputs: Map[DxName, JsValue]): Map[DxName, Value] = {
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
