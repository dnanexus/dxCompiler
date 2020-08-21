package dx.exec

import dx.api.{DxApi, DxFile}
import dx.core.languages.wdl.{WdlExecutableLink, TypeSerialization, WdlDxLinkSerde}
import spray.json._
import wdlTools.eval.WdlValues
import wdlTools.types.WdlTypes

case class WfFragInput(blockPath: Vector[Int],
                       env: Map[String, (WdlTypes.T, WdlValues.V)],
                       execLinkInfo: Map[String, WdlExecutableLink])

case class WfFragInputOutput(typeAliases: Map[String, WdlTypes.T],
                             wdlVarLinksConverter: WdlDxLinkSerde,
                             dxApi: DxApi) {

  private def loadWorkflowMetaInfo(
      metaInfo: Map[String, JsValue]
  ): (Map[String, WdlExecutableLink], Vector[Int], Map[String, WdlTypes.T]) = {
    // meta information used for running workflow fragments
    val execLinkInfo: Map[String, WdlExecutableLink] = metaInfo.get("execLinkInfo") match {
      case None => Map.empty
      case Some(JsObject(fields)) =>
        fields.map {
          case (key, ali) =>
            key -> WdlExecutableLink.readJson(dxApi, ali, typeAliases)
        }
      case other => throw new Exception(s"Bad value ${other}")
    }
    val blockPath: Vector[Int] = metaInfo.get("blockPath") match {
      case None => Vector.empty
      case Some(JsArray(arr)) =>
        arr.map {
          case JsNumber(n) => n.toInt
          case _           => throw new Exception("Bad value ${arr}")
        }
      case other => throw new Exception(s"Bad value ${other}")
    }
    val fqnDictTypes: Map[String, WdlTypes.T] = metaInfo.get("fqnDictTypes") match {
      case Some(JsObject(fields)) =>
        fields.map {
          case (key, JsString(value)) =>
            // Transform back to a fully qualified name with dots
            val orgKeyName = WdlDxLinkSerde.decodeDots(key)
            val wdlType = TypeSerialization(typeAliases).fromString(value)
            orgKeyName -> wdlType
          case other => throw new Exception(s"Bad value ${other}")
        }
      case other => throw new Exception(s"Bad value ${other}")
    }

    (execLinkInfo, blockPath, fqnDictTypes)
  }

  // 1. Convert the inputs to WDL values
  // 2. Setup an environment to evaluate the sub-block. This should
  //    look to the WDL code as if all previous code had been evaluated.
  def loadInputs(inputs: JsValue, metaInfo: JsValue): WfFragInput = {
    val regularFields: Map[String, JsValue] = inputs.asJsObject.fields
      .filter { case (fieldName, _) => !fieldName.endsWith(WdlDxLinkSerde.FlatFilesSuffix) }

    // Extract the meta information needed to setup the closure for the subblock
    val (execLinkInfo, blockPath, fqnDictTypes) = loadWorkflowMetaInfo(metaInfo.asJsObject.fields)

    // What remains are inputs from other stages. Convert from JSON to WDL values
    val env: Map[String, (WdlTypes.T, WdlValues.V)] = regularFields.map {
      case (name, jsValue) =>
        val fqn = WdlDxLinkSerde.decodeDots(name)
        val wdlType = fqnDictTypes.get(fqn) match {
          case None =>
            throw new Exception(s"Did not find variable ${fqn} (${name}) in the block environment")
          case Some(x) => x
        }
        val value = wdlVarLinksConverter.unpackJobInput(fqn, wdlType, jsValue)
        fqn -> (wdlType, value)
    }

    WfFragInput(blockPath, env, execLinkInfo)
  }

  // find all the dx:files that are referenced from the inputs
  def findRefDxFiles(inputs: JsValue, metaInfo: JsValue): Vector[DxFile] = {
    val regularFields: Map[String, JsValue] = inputs.asJsObject.fields
      .filter { case (fieldName, _) => !fieldName.endsWith(WdlDxLinkSerde.FlatFilesSuffix) }

    val (_, _, fqnDictTypes) = loadWorkflowMetaInfo(metaInfo.asJsObject.fields)

    // Convert from JSON to WDL values
    regularFields
      .map {
        case (name, jsValue) =>
          val fqn = WdlDxLinkSerde.decodeDots(name)
          if (!fqnDictTypes.contains(fqn)) {
            throw new Exception(
                s"Did not find variable ${fqn} (${name}) in the block environment"
            )
          }
          dxApi.findFiles(jsValue)
      }
      .toVector
      .flatten
  }
}
