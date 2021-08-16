package dx

import dx.core.ir.Type.{TArray, TFile, THash, TOptional, TSchema}
import dx.core.ir.{Type, Value}
import dx.core.ir.Value.{VArray, VFile, VHash, VNull, VString}
import dx.util.{AddressableFileNode, FileSourceResolver, LocalFileSource}

import java.nio.file.Files

package object executor {
  // TODO: it would be nice to extract dx:// links from VString values - this will
  //  happen in the case where the container is a dx file and being passed in as
  //  an input parameter - so that they could be downloaded using dxda. However,
  //  this would also require some way for the downloaded image tarball to be
  //  discovered and loaded. For now, we rely on DockerUtils to download the image
  //  (via DxFileSource, which uses the dx API to download the file).
  def extractFiles(
      v: Value,
      fileResolver: FileSourceResolver
  ): Vector[AddressableFileNode] = {
    def inner(innerValue: Value): Vector[AddressableFileNode] = {
      innerValue match {
        case VFile(s)      => Vector(fileResolver.resolve(s))
        case VArray(value) => value.flatMap(inner)
        case VHash(m)      => m.values.flatMap(inner).toVector
        case _             => Vector.empty
      }
    }
    inner(v)
  }

  def extractOutputFiles(
      name: String,
      v: Value,
      t: Type,
      fileResolver: FileSourceResolver,
      ignoreNonLocalFiles: Boolean = false
  ): Vector[AddressableFileNode] = {
    def getFileNode(varName: String,
                    fs: AddressableFileNode,
                    optional: Boolean): Vector[AddressableFileNode] = {
      fs match {
        case localFs: LocalFileSource if optional && !Files.exists(localFs.canonicalPath) =>
          // ignore optional, non-existent files
          Vector.empty
        case localFs: LocalFileSource if !Files.exists(localFs.canonicalPath) =>
          throw new Exception(
              s"required output file ${varName} does not exist at ${localFs.canonicalPath}"
          )
        case localFs: LocalFileSource     => Vector(localFs)
        case other if ignoreNonLocalFiles => Vector.empty
        case other =>
          throw new RuntimeException(s"${varName} specifies non-local file ${other}")
      }
    }

    def inner(innerName: String,
              innerValue: Value,
              innerType: Type,
              optional: Boolean): Vector[AddressableFileNode] = {
      (innerType, innerValue) match {
        case (TOptional(_), VNull) =>
          Vector.empty
        case (TOptional(optType), _) =>
          inner(innerName, innerValue, optType, optional = true)
        case (TFile, VFile(path)) =>
          getFileNode(innerName, fileResolver.resolve(path), optional)
        case (TFile, VString(path)) =>
          getFileNode(innerName, fileResolver.resolve(path), optional)
        case (TArray(_, nonEmpty), VArray(array)) if nonEmpty && array.isEmpty =>
          throw new Exception(s"Non-empty array ${name} has empty value")
        case (TArray(elementType, _), VArray(array)) =>
          array.zipWithIndex.flatMap {
            case (element, index) =>
              inner(s"${innerName}[${index}]", element, elementType, optional = false)
          }
        case (TSchema(name, memberTypes), VHash(members)) =>
          memberTypes.flatMap {
            case (key, t) =>
              (t, members.get(key)) match {
                case (TOptional(_), None) =>
                  Vector.empty
                case (_, None) =>
                  throw new Exception(s"missing non-optional member ${key} of struct ${name}")
                case (_, Some(v)) =>
                  inner(s"${name}.${key}", v, t, optional = false)
              }
          }.toVector
        case (THash, VHash(members)) =>
          members.flatMap {
            case (key, value) =>
              val files = extractFiles(value, fileResolver)
              files.flatMap(fs => getFileNode(s"${innerName}.${key}", fs, optional = true))
          }.toVector
        case _ =>
          Vector.empty
      }
    }

    inner(name, v, t, optional = false)
  }

  def delocalizeOutputFiles(localizedOutputs: Map[String, (Type, Value)],
                            resolveUri: String => Option[String],
                            ignoreMissingUris: Boolean = false): Map[String, (Type, Value)] = {
    def pathTranslator(v: Value, t: Option[Type], optional: Boolean): Option[Value] = {
      val uri = (t, v) match {
        case (_, VFile(uri))             => Some(uri)
        case (Some(TFile), VString(uri)) => Some(uri)
        case _                           => None
      }
      uri.flatMap { u =>
        resolveUri(u) match {
          case Some(uri)                 => Some(VFile(uri))
          case None if ignoreMissingUris => None
          case None if optional          => Some(VNull)
          case None =>
            throw new Exception(s"Did not delocalize file ${u}")
        }
      }
    }
    localizedOutputs.view.mapValues {
      case (t, v) => (t, Value.transform(v, Some(t), pathTranslator))
    }.toMap
  }
}
