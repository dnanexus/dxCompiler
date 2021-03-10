package dx.core.languages.wdl

import dx.api.DxPath

import java.nio.file.Path
import dx.core.ir.{Type, TypeSerde, Value}
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.util.{Bindings, Enum, FileNode, FileSourceResolver, Logger, StringFileNode}
import wdlTools.eval.{Coercion, EvalUtils, VBindings}
import wdlTools.eval.WdlValues._
import wdlTools.syntax.{
  NoSuchParserException,
  Parsers,
  SourceLocation,
  SyntaxException,
  WdlParser,
  WdlVersion,
  AbstractSyntax => AST
}
import wdlTools.types.{
  Section,
  TypeCheckingRegime,
  TypeException,
  TypeInfer,
  TypeUtils,
  TypedAbstractSyntax => TAT
}
import wdlTools.types.WdlTypes._

import scala.collection.immutable.{SeqMap, TreeSeqMap}

/**
  * The kind of block input variable being referenced. A Computed
  * input is one that is computed from other inputs - currently
  * this is only used for the scatter variable.
  */
object InputKind extends Enum {
  type InputKind = Value
  val Required, Computed, Optional = Value
}

/**
  * A reference to a block input variable.
  * @param identifierParts the parts of the identifier name
  * @param fieldName the fieldName within the call/object
  *                  referred to by the identifier, if any
  * @param wdlType the variable's type
  * @param kind the variable's kind
  */
private case class WdlInputRef(identifierParts: Vector[String],
                               fieldName: Option[String],
                               wdlType: T,
                               kind: InputKind.InputKind) {
  lazy val identifier: String = identifierParts.mkString(".")

  lazy val fullyQualifiedName: String = {
    (identifierParts ++ fieldName.toVector).mkString(".")
  }

  /**
    * Returns an Iterator over all the nested identifiers
    * represented by this reference.
    */
  def nameIter: Iterator[String] = {
    (identifierParts ++ fieldName.toVector).reverse.tails.collect {
      case parts if parts.nonEmpty => parts.reverse.mkString(".")
    }
  }
}

object WdlUtils {

  val locPlaceholder: SourceLocation = SourceLocation.empty

  def parseSource(sourceCode: FileNode,
                  parser: WdlParser,
                  logger: Logger = Logger.get): AST.Document = {
    try {
      parser.parseDocument(sourceCode)
    } catch {
      case se: SyntaxException =>
        logger.error(s"WDL code is syntactically invalid -----\n${sourceCode.readString}")
        throw se
    }
  }

  def parseAndCheckSource(
      sourceCode: FileNode,
      parser: WdlParser,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      regime: TypeCheckingRegime.TypeCheckingRegime = TypeCheckingRegime.Moderate,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, T_Struct]) = {
    val doc = parseSource(sourceCode, parser)
    try {
      val (tDoc, ctx) =
        TypeInfer(regime, fileResolver = fileResolver, logger = logger)
          .apply(doc)
      (tDoc, ctx.aliases)
    } catch {
      case te: TypeException =>
        logger.error(
            s"WDL code is syntactically valid BUT it fails type-checking -----\n${sourceCode.readString}"
        )
        throw te
    }
  }

  def parseAndCheckSourceNode(
      node: FileNode,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      regime: TypeCheckingRegime.TypeCheckingRegime = TypeCheckingRegime.Moderate,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, T_Struct]) = {
    val parser =
      try {
        Parsers(followImports = true, fileResolver = fileResolver, logger = logger).getParser(node)
      } catch {
        case nspe: NoSuchParserException =>
          logger.error(
              s"Source code does not appear to be any supported version of WDL -----\n${node.readString}"
          )
          throw nspe
      }
    parseAndCheckSource(node, parser, fileResolver, regime, logger)
  }

  /**
    * Parses a top-level WDL file and all its imports.
    * @param path the path to the WDL file
    * @param fileResolver FileSourceResolver
    * @param regime TypeCheckingRegime
    * @param logger Logger
    * @return (document, aliases), where aliases is a mapping of all the (fully-qualified)
    *         alias names to values. Aliases include Structs defined in any file (which are
    *         *not* namespaced) and all aliases defined in import statements of all documents
    *         (which *are* namespaced).
    */
  def parseAndCheckSourceFile(
      path: Path,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      regime: TypeCheckingRegime.TypeCheckingRegime = TypeCheckingRegime.Moderate,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, T_Struct]) = {
    parseAndCheckSourceNode(fileResolver.fromPath(path), fileResolver, regime, logger)
  }

  def parseAndCheckSourceString(
      sourceCodeStr: String,
      name: String,
      fileResolver: FileSourceResolver = FileSourceResolver.get,
      regime: TypeCheckingRegime.TypeCheckingRegime = TypeCheckingRegime.Moderate,
      logger: Logger = Logger.get
  ): (TAT.Document, Bindings[String, T_Struct]) = {
    parseAndCheckSourceNode(StringFileNode(sourceCodeStr, name), fileResolver, regime, logger)
  }

  def parseExpr(exprStr: String,
                wdlVersion: WdlVersion,
                parser: WdlParser,
                docSource: FileNode,
                bindings: Bindings[String, T],
                section: Section.Section = Section.Other,
                fileResolver: FileSourceResolver = FileSourceResolver.get,
                regime: TypeCheckingRegime.TypeCheckingRegime = TypeCheckingRegime.Moderate,
                logger: Logger = Logger.get): TAT.Expr = {
    val expr = parser.parseExpr(exprStr)
    try {
      TypeInfer(regime, fileResolver = fileResolver, logger = logger)
        .applyExpression(expr, wdlVersion, regime, docSource, bindings, section)
    } catch {
      case te: TypeException =>
        logger.error(
            s"Expression code is syntactically valid BUT it fails type-checking: ${exprStr}"
        )
        throw te
    }
  }

  def getUnqualifiedName(name: String): String = {
    if (name contains ".") {
      name.split("\\.").last
    } else {
      name
    }
  }

  // create a wdl-value of a specific type.
  def getDefaultValueOfType(wdlType: T, loc: SourceLocation = WdlUtils.locPlaceholder): TAT.Expr = {
    wdlType match {
      case T_Boolean => TAT.ValueBoolean(value = true, wdlType, loc)
      case T_Int     => TAT.ValueInt(0, wdlType, loc)
      case T_Float   => TAT.ValueFloat(0.0, wdlType, loc)
      case T_String  => TAT.ValueString("", wdlType, loc)
      case T_File    => TAT.ValueString("placeholder.txt", wdlType, loc)

      // We could convert an optional to a null value, but that causes
      // problems for the pretty printer.
      // WdlValues.V_OptionalValue(wdlType, None)
      case T_Optional(t) => getDefaultValueOfType(t)

      // The WdlValues.V_Map type HAS to appear before the array types, because
      // otherwise it is coerced into an array. The map has to
      // contain at least one key-value pair, otherwise you get a type error.
      case T_Map(keyType, valueType) =>
        val k = getDefaultValueOfType(keyType)
        val v = getDefaultValueOfType(valueType)
        TAT.ExprMap(TreeSeqMap(k -> v), wdlType, loc)

      // an empty array
      case T_Array(_, false) =>
        TAT.ExprArray(Vector.empty, wdlType, loc)

      // Non empty array
      case T_Array(t, true) =>
        TAT.ExprArray(Vector(getDefaultValueOfType(t)), wdlType, loc)

      case T_Pair(lType, rType) =>
        TAT.ExprPair(getDefaultValueOfType(lType), getDefaultValueOfType(rType), wdlType, loc)

      case T_Struct(_, typeMap) =>
        val members = typeMap.map {
          case (fieldName, t) =>
            val key: TAT.Expr = TAT.ValueString(fieldName, T_String, loc)
            key -> getDefaultValueOfType(t)
        }
        TAT.ExprObject(members, wdlType, loc)

      case T_Object =>
        TAT.ExprObject(SeqMap.empty, wdlType, SourceLocation.empty)

      case _ => throw new Exception(s"Unhandled type ${wdlType}")
    }
  }

  // Functions to convert between WDL and IR types and values.
  // WDL has two types that IR does not: Pair and Map. These
  // are represented as objects with specific keys ('left' and
  // 'right' for Pair, 'keys' and 'values' for Map). We define
  // a special TSchema for each of these types.

  /**
    * Name prefix for Pair-type schemas.
    */
  val PairSchemaPrefix = "Pair___"

  /**
    * Pair.left key that can be used in WDL input files.
    */
  val PairLeftKey = "left"

  /**
    * Pair.right key that can be used in WDL input files.
    */
  val PairRightKey = "right"

  def createPairSchema(left: Type, right: Type): TSchema = {
    val name = s"${PairSchemaPrefix}(${TypeSerde.toString(left)}, ${TypeSerde.toString(right)})"
    TSchema(name, TreeSeqMap(PairLeftKey -> left, PairRightKey -> right))
  }

  def isPairSchema(t: TSchema): Boolean = {
    t.name.startsWith(PairSchemaPrefix) && t.fields.size == 2 && t.fields.keySet == Set(
        PairLeftKey,
        PairRightKey
    )
  }

  def isPairValue(fields: Map[String, _]): Boolean = {
    fields.size == 2 && fields.keySet == Set(PairLeftKey, PairRightKey)
  }

  /**
    * Name prefix for Map-type schemas.
    */
  val MapSchemaPrefix = "Map___"

  /**
    * Map.keys key that can be used in WDL input files.
    */
  val MapKeysKey = "keys"

  /**
    * Map.values key that can be used in WDL  input files.
    */
  val MapValuesKey = "values"

  def createMapSchema(keyType: Type, valueType: Type): TSchema = {
    val name =
      s"${MapSchemaPrefix}[${TypeSerde.toString(keyType)}, ${TypeSerde.toString(valueType)}]"
    TSchema(name, TreeSeqMap(MapKeysKey -> TArray(keyType), MapValuesKey -> TArray(valueType)))
  }

  def isMapSchema(t: TSchema): Boolean = {
    t.name.startsWith(MapSchemaPrefix) && t.fields.size == 2 && t.fields.keySet == Set(
        MapKeysKey,
        MapValuesKey
    )
  }

  def isMapValue(fields: Map[String, _]): Boolean = {
    fields.size == 2 && fields.keySet == Set(MapKeysKey, MapValuesKey)
  }

  def mapValueToMap(fields: Map[String, Value]): Map[Value, Value] = {
    val keys = fields.get(MapKeysKey) match {
      case Some(VArray(array)) => array
      case _                   => throw new Exception(s"invalid map value ${fields}")
    }
    val values = fields.get(MapValuesKey) match {
      case Some(VArray(array)) => array
      case _                   => throw new Exception(s"invalid map value ${fields}")
    }
    keys.zip(values).toMap
  }

  def toIRSchema(wdlStruct: T_Struct): TSchema = {
    TSchema(wdlStruct.name, wdlStruct.members.map {
      case (key, value) => key -> toIRType(value)
    })
  }
  def toIRType(wdlType: T): Type = {
    wdlType match {
      case T_Boolean     => TBoolean
      case T_Int         => TInt
      case T_Float       => TFloat
      case T_String      => TString
      case T_File        => TFile
      case T_Directory   => TDirectory
      case T_Object      => THash
      case T_Optional(t) => TOptional(toIRType(t))
      case T_Array(t, nonEmpty) =>
        TArray(toIRType(t), nonEmpty)
      case struct: T_Struct =>
        toIRSchema(struct)
      case T_Pair(leftType, rightType) =>
        createPairSchema(toIRType(leftType), toIRType(rightType))
      case T_Map(keyType, valueType) =>
        createMapSchema(toIRType(keyType), toIRType(valueType))
      case _ =>
        throw new Exception(s"Cannot convert WDL type ${wdlType} to IR")
    }
  }

  def toIRTypeMap(wdlTypes: Map[String, T]): Map[String, Type] = {
    wdlTypes.map {
      case (name, t) => name -> toIRType(t)
    }
  }

  def toIRSchemaMap(wdlTypes: Map[String, T_Struct]): Map[String, TSchema] = {
    wdlTypes.map {
      case (name, t) => name -> toIRSchema(t)
    }
  }

  def fromIRType(irType: Type, typeAliases: Map[String, T] = Map.empty): T = {
    def inner(innerType: Type): T = {
      innerType match {
        case TBoolean     => T_Boolean
        case TInt         => T_Int
        case TFloat       => T_Float
        case TString      => T_String
        case TFile        => T_File
        case TDirectory   => T_Directory
        case THash        => T_Object
        case TOptional(t) => T_Optional(inner(t))
        case TArray(t, nonEmpty) =>
          T_Array(inner(t), nonEmpty = nonEmpty)
        case TSchema(name, _) if typeAliases.contains(name) =>
          typeAliases(name)
        case pairSchema: TSchema if isPairSchema(pairSchema) =>
          T_Pair(inner(pairSchema.fields(PairLeftKey)), inner(pairSchema.fields(PairRightKey)))
        case mapSchema: TSchema if isMapSchema(mapSchema) =>
          val keyType = inner(mapSchema.fields(MapKeysKey)) match {
            case T_Array(t, _) => t
            case other         => throw new Exception(s"invalid Map schema key type ${other}")
          }
          val valueType = inner(mapSchema.fields(MapValuesKey)) match {
            case T_Array(t, _) => t
            case other         => throw new Exception(s"invalid Map schema value type ${other}")
          }
          T_Map(keyType, valueType)
        case TSchema(name, _) =>
          throw new Exception(s"Unknown type ${name}")
        case _ =>
          throw new Exception(s"Cannot convert IR type ${innerType} to WDL")
      }
    }
    inner(irType)
  }

  def toIRValue(wdlValue: V): Value = {
    wdlValue match {
      case V_Null            => VNull
      case V_Boolean(b)      => VBoolean(value = b)
      case V_Int(i)          => VInt(i)
      case V_Float(f)        => VFloat(f)
      case V_String(s)       => VString(s)
      case V_File(path)      => VFile(path)
      case V_Directory(path) => VFolder(path)
      case V_Array(array) =>
        VArray(array.map(v => toIRValue(v)))
      case V_Pair(left, right) =>
        // encode this as a hash with 'left' and 'right' keys
        VHash(
            TreeSeqMap(
                PairLeftKey -> toIRValue(left),
                PairRightKey -> toIRValue(right)
            )
        )
      case V_Map(members) =>
        // encode this as a hash with 'keys' and 'values' keys
        val (keys, values) = members.map {
          case (k, v) => (toIRValue(k), toIRValue(v))
        }.unzip
        VHash(
            TreeSeqMap(
                MapKeysKey -> VArray(keys.toVector),
                MapValuesKey -> VArray(values.toVector)
            )
        )
      case V_Object(fields) =>
        VHash(fields.map {
          case (key, value) => key -> toIRValue(value)
        })
      case V_Struct(_, fields) =>
        VHash(fields.map {
          case (key, value) => key -> toIRValue(value)
        })
      case _ =>
        throw new Exception(s"Invalid WDL value ${wdlValue})")
    }
  }

  def toIRValue(wdlValue: V, wdlType: T): Value = {
    (wdlType, wdlValue) match {
      case (_, V_Null)                      => VNull
      case (T_Boolean, V_Boolean(b))        => VBoolean(value = b)
      case (T_Int, V_Int(i))                => VInt(i)
      case (T_Float, V_Float(f))            => VFloat(f)
      case (T_String, V_String(s))          => VString(s)
      case (T_File, V_String(path))         => VFile(path)
      case (T_File, V_File(path))           => VFile(path)
      case (T_Directory, V_String(path))    => VFolder(path)
      case (T_Directory, V_Directory(path)) => VFolder(path)
      case (T_Object, o: V_Object)          => toIRValue(o)
      case (T_Optional(t), V_Optional(v))   => toIRValue(v, t)
      case (T_Optional(t), v)               => toIRValue(v, t)
      case (t, V_Optional(v))               => toIRValue(v, t)
      case (T_Array(_, true), V_Array(array)) if array.isEmpty =>
        throw new Exception(
            s"Empty array with non-empty (+) quantifier"
        )
      case (T_Array(t, _), V_Array(array)) =>
        VArray(array.map(v => toIRValue(v, t)))
      case (T_Pair(leftType, rightType), V_Pair(leftValue, rightValue)) =>
        // encode this as a hash with left and right keys
        VHash(
            TreeSeqMap(
                PairLeftKey -> toIRValue(leftValue, leftType),
                PairRightKey -> toIRValue(rightValue, rightType)
            )
        )
      case (T_Map(keyType, valueType), V_Map(members)) =>
        // encode this as a hash with 'keys' and 'values' keys
        val (keys, values) = members.map {
          case (k, v) => (toIRValue(k, keyType), toIRValue(v, valueType))
        }.unzip
        VHash(
            TreeSeqMap(
                MapKeysKey -> VArray(keys.toVector),
                MapValuesKey -> VArray(values.toVector)
            )
        )
      case (T_Struct(name, memberTypes), V_Struct(vName, memberValues)) if name == vName =>
        structToIRValue(name, memberValues, memberTypes)
      case (T_Struct(name, memberTypes), V_Object(memberValues)) =>
        structToIRValue(name, memberValues, memberTypes)
      case _ =>
        throw new Exception(s"Invalid (type, value) combination (${wdlType}, ${wdlValue})")
    }
  }

  private def structToIRValue(name: String,
                              fieldValues: SeqMap[String, V],
                              fieldTypes: SeqMap[String, T]): VHash = {
    // ensure 1) members keys are a subset of memberTypes keys, 2) members
    // values are convertable to the corresponding types, and 3) any keys
    // in memberTypes that do not appear in members are optional
    val keys1 = fieldValues.keySet
    val keys2 = fieldTypes.keySet
    val extra = keys2.diff(keys1)
    if (extra.nonEmpty) {
      throw new Exception(
          s"struct ${name} value has members that do not appear in the struct definition: ${extra}"
      )
    }
    val missingNonOptional = keys1.diff(keys2).map(key => key -> fieldTypes(key)).filterNot {
      case (_, T_Optional(_)) => false
      case _                  => true
    }
    if (missingNonOptional.nonEmpty) {
      throw new Exception(
          s"struct ${name} value is missing non-optional members ${missingNonOptional}"
      )
    }
    VHash(fieldTypes.collect {
      case (name, t) if fieldValues.contains(name) =>
        name -> toIRValue(fieldValues(name), t)
    })
  }

  def toIR(wdl: Map[String, (T, V)]): Map[String, (Type, Value)] = {
    wdl.map {
      case (name, (wdlType, wdlValue)) =>
        val irType = toIRType(wdlType)
        val irValue = toIRValue(wdlValue, wdlType)
        name -> (irType, irValue)
    }
  }

  def fromIRValue(value: Value, name: Option[String]): V = {
    value match {
      case VNull       => V_Null
      case VBoolean(b) => V_Boolean(b)
      case VInt(i)     => V_Int(i)
      case VFloat(f)   => V_Float(f)
      case VString(s)  => V_String(s)
      case f: VFile    => V_File(f.uri)
      case f: VFolder  => V_Directory(f.uri)
      case VArray(array) =>
        V_Array(array.zipWithIndex.map {
          case (v, i) => fromIRValue(v, name.map(n => s"${n}[${i}]"))
        })
      case VHash(fields) if isPairValue(fields) =>
        V_Pair(
            fromIRValue(fields(PairLeftKey), name.map(n => s"${n}.${PairLeftKey}")),
            fromIRValue(fields(PairRightKey), name.map(n => s"${n}.${PairRightKey}"))
        )
      case VHash(fields) if isMapValue(fields) =>
        val keys = fromIRValue(fields(MapKeysKey), name.map(n => s"${n}[${MapKeysKey}]"))
        val values =
          fromIRValue(fields(MapValuesKey), name.map(n => s"${n}[${MapValuesKey}]"))
        (keys, values) match {
          case (V_Array(keyArray), V_Array(valueArray)) =>
            V_Map(keyArray.zip(valueArray).to(TreeSeqMap))
          case other =>
            throw new Exception(s"invalid map value ${other}")
        }
      case VHash(fields) =>
        V_Object(fields.map {
          case (key, value) => key -> fromIRValue(value, name.map(n => s"${n}[${key}]"))
        })
      case _ =>
        throw new Exception(
            s"cannot convert ${name.getOrElse("IR")} value ${value} to WDL value"
        )
    }
  }

  def fromIRValue(value: Value, wdlType: T, name: String): V = {
    (wdlType, value) match {
      case (T_Optional(_), VNull)       => V_Null
      case (T_Boolean, VBoolean(b))     => V_Boolean(value = b)
      case (T_Int, VInt(i))             => V_Int(i)
      case (T_Float, VFloat(f))         => V_Float(f)
      case (T_Float, VInt(i))           => V_Float(i.toDouble)
      case (T_String, VString(s))       => V_String(s)
      case (T_File, VString(path))      => V_File(path)
      case (T_File, f: VFile)           => V_File(f.uri)
      case (T_Directory, VString(path)) => V_Directory(path)
      case (T_Directory, f: VFolder)    => V_Directory(f.uri)
      case (T_Object, o: VHash)         => fromIRValue(o, Some(name))
      case (T_Optional(t), v)           => V_Optional(fromIRValue(v, t, name))
      case (T_Array(_, true), VArray(array)) if array.isEmpty =>
        throw new Exception(
            s"Empty array with non-empty (+) quantifier"
        )
      case (T_Array(t, _), VArray(array)) =>
        V_Array(array.zipWithIndex.map {
          case (v, i) => fromIRValue(v, t, s"${name}[${i}]")
        })
      case (T_Pair(leftType, rightType), VHash(fields)) if isPairValue(fields) =>
        V_Pair(
            fromIRValue(fields(PairLeftKey), leftType, s"${name}.${PairLeftKey}"),
            fromIRValue(fields(PairRightKey), rightType, s"${name}.${PairRightKey}")
        )
      case (T_Map(keyType, valueType), VHash(fields)) if isMapValue(fields) =>
        // keyType and valueType will be the map element types, but the keys
        // and values are encoded as arrays so we need to wrap the types in T_Array
        val keys = fromIRValue(fields(MapKeysKey), T_Array(keyType), s"${name}[${MapKeysKey}]")
        val values =
          fromIRValue(fields(MapValuesKey), T_Array(valueType), s"${name}[${MapValuesKey}]")
        (keys, values) match {
          case (V_Array(keyArray), V_Array(valueArray)) =>
            V_Map(keyArray.zip(valueArray).to(TreeSeqMap))
          case other =>
            throw new Exception(s"invalid map value ${other}")
        }
      case (T_Struct(structName, memberTypes), VHash(members)) =>
        // ensure 1) members keys are a subset of memberTypes keys, 2) members
        // values are convertable to the corresponding types, and 3) any keys
        // in memberTypes that do not appear in members are optional
        val keys1 = members.keySet
        val keys2 = memberTypes.keySet
        val extra = keys2.diff(keys1)
        if (extra.nonEmpty) {
          throw new Exception(
              s"struct ${structName} value has members that do not appear in the struct definition: ${extra}"
          )
        }
        val missingNonOptional = keys1.diff(keys2).map(key => key -> memberTypes(key)).filterNot {
          case (_, T_Optional(_)) => false
          case _                  => true
        }
        if (missingNonOptional.nonEmpty) {
          throw new Exception(
              s"struct ${structName} value is missing non-optional members ${missingNonOptional}"
          )
        }
        V_Object(members.map {
          case (key, value) => key -> fromIRValue(value, memberTypes(key), s"${name}[${key}]")
        })
      case _ =>
        throw new Exception(
            s"Cannot convert ${name} (${wdlType}, ${value}) to WDL value"
        )
    }
  }

  def fromIR(ir: Map[String, (Type, Value)],
             typeAliases: Map[String, T] = Map.empty): Map[String, (T, V)] = {
    ir.map {
      case (name, (t, v)) =>
        val wdlType = fromIRType(t, typeAliases)
        val wdlValue = fromIRValue(v, wdlType, name)
        name -> (wdlType, wdlValue)
    }
  }

  private def ensureUniformType(exprs: Iterable[TAT.Expr]): T = {
    exprs.headOption.map(_.wdlType) match {
      case Some(t) if exprs.tail.exists(_.wdlType != t) =>
        throw new Exception(s"${exprs} contains non-homogeneous values")
      case Some(t) => t
      case None    => T_Any
    }
  }

  def irValueToExpr(value: Value): TAT.Expr = {
    val loc = SourceLocation.empty
    value match {
      case VNull       => TAT.ValueNull(T_Any, loc)
      case VBoolean(b) => TAT.ValueBoolean(b, T_Boolean, loc)
      case VInt(i)     => TAT.ValueInt(i, T_Int, loc)
      case VFloat(f)   => TAT.ValueFloat(f, T_Float, loc)
      case VString(s)  => TAT.ValueString(s, T_String, loc)
      case f: VFile    => TAT.ValueFile(f.uri, T_File, loc)
      case f: VFolder  => TAT.ValueDirectory(f.uri, T_Directory, loc)
      case VArray(array) =>
        val a = array.map(irValueToExpr)
        val t = ensureUniformType(a)
        TAT.ExprArray(a, t, loc)
      case VHash(fields) if isPairValue(fields) =>
        val left = irValueToExpr(fields(PairLeftKey))
        val right = irValueToExpr(fields(PairRightKey))
        TAT.ExprPair(left, right, T_Pair(left.wdlType, right.wdlType), loc)
      case VHash(fields) if isMapValue(fields) =>
        val keys = irValueToExpr(fields(MapKeysKey))
        val values = irValueToExpr(fields(MapValuesKey))
        (keys, values) match {
          case (TAT.ExprArray(keyArray, keyType, _), TAT.ExprArray(valueArray, valueType, _)) =>
            TAT.ExprMap(keyArray.zip(valueArray).to(TreeSeqMap), T_Map(keyType, valueType), loc)
          case other =>
            throw new Exception(s"invalid map value ${other}")
        }
      case VHash(members) =>
        val m: SeqMap[TAT.Expr, TAT.Expr] = members
          .map {
            case (key, value) => TAT.ValueString(key, T_String, loc) -> irValueToExpr(value)
          }
          .to(TreeSeqMap)
        TAT.ExprObject(m, T_Object, loc)
      case _ =>
        throw new Exception(s"cannot convert IR value ${value} to WDL")
    }
  }

  /**
    * A trivial expression has no operators, it is either(1) a constant,
    * (2) a single identifier, or (3) an access to a call field.
    * For example, `5`, `['a', 'b', 'c']`, and `true` are trivial.
    * 'x + y'  is not.
    * @param expr expression
    * @return
    */
  def isTrivialExpression(expr: TAT.Expr): Boolean = {
    expr match {
      case expr if TypeUtils.isPrimitiveValue(expr) => true
      case _: TAT.ExprIdentifier                    => true

      // A collection of constants
      case TAT.ExprPair(l, r, _, _)   => Vector(l, r).forall(TypeUtils.isPrimitiveValue)
      case TAT.ExprArray(value, _, _) => value.forall(TypeUtils.isPrimitiveValue)
      case TAT.ExprMap(value, _, _) =>
        value.forall {
          case (k, v) => TypeUtils.isPrimitiveValue(k) && TypeUtils.isPrimitiveValue(v)
        }
      case TAT.ExprObject(value, _, _) => value.values.forall(TypeUtils.isPrimitiveValue)

      // Access a field in a call or a struct
      //   Int z = eliminateDuplicate.fields
      case TAT.ExprGetName(_: TAT.ExprIdentifier, _, _, _) => true

      case _ => false
    }
  }

  def isDxFile(file: V_File): Boolean = {
    file.value.startsWith(DxPath.DxUriPrefix)
  }

  /**
    * Deep search for all calls in WorkflowElements.
    * @param elements WorkflowElements
    * @return
    */
  def deepFindCalls(elements: Vector[TAT.WorkflowElement]): Vector[TAT.Call] = {
    elements.foldLeft(Vector.empty[TAT.Call]) {
      case (accu, call: TAT.Call) =>
        accu :+ call
      case (accu, ssc: TAT.Scatter) =>
        accu ++ deepFindCalls(ssc.body)
      case (accu, ifStmt: TAT.Conditional) =>
        accu ++ deepFindCalls(ifStmt.body)
      case (accu, _) =>
        accu
    }
  }

  /**
    * Returns the variables used in an expression. Nested references do
    * not include members, except when the identifier refers to a call or
    * `withField` is true.
    * @param expr the expression
    * @param withField whether to append the field name of a terminal
    *                  GetName expression to the LHS (requires the LHS
    *                  to be a Pair, object, or struct; always true if
    *                  the LHS is a call).
    * @return Vector of (name, type, optional) tuples
    * @example
    * expression   inputs
    *   1 + 2        Vector.empty
    *   x + y        Vector(('x', None), ('y', None))
    *   foo.y        Vector(('foo', Some('y')))    [foo is a call]
    *   bar.left     Vector(('bar', None))         [bar is a Pair, withField = false]
    *   bar.left     Vector(('bar', Some('left'))) [bar is a Pair, withField = true]
    */
  private def getExpressionInputs(expr: TAT.Expr, withField: Boolean): Vector[WdlInputRef] = {
    def inner(
        innerExpr: TAT.Expr
    ): Vector[WdlInputRef] = {
      innerExpr match {
        case _: TAT.ValueNull      => Vector.empty
        case _: TAT.ValueNone      => Vector.empty
        case _: TAT.ValueBoolean   => Vector.empty
        case _: TAT.ValueInt       => Vector.empty
        case _: TAT.ValueFloat     => Vector.empty
        case _: TAT.ValueString    => Vector.empty
        case _: TAT.ValueFile      => Vector.empty
        case _: TAT.ValueDirectory => Vector.empty
        case TAT.ExprIdentifier(id, wdlType, _) =>
          val kind = if (TypeUtils.isOptional(wdlType)) {
            InputKind.Optional
          } else {
            InputKind.Required
          }
          Vector(WdlInputRef(Vector(id), None, wdlType, kind))
        case TAT.ExprCompoundString(valArr, _, _) =>
          valArr.flatMap(elem => inner(elem))
        case TAT.ExprPair(l, r, _, _) =>
          inner(l) ++ inner(r)
        case TAT.ExprArray(arrVal, _, _) =>
          arrVal.flatMap(elem => inner(elem))
        case TAT.ExprMap(valMap, _, _) =>
          valMap
            .map { case (k, v) => inner(k) ++ inner(v) }
            .toVector
            .flatten
        case TAT.ExprObject(fields, _, _) =>
          fields
            .map { case (_, v) => inner(v) }
            .toVector
            .flatten
        case TAT.ExprPlaceholderCondition(t: TAT.Expr, f: TAT.Expr, value: TAT.Expr, _, _) =>
          inner(t) ++ inner(f) ++ inner(value)
        case TAT.ExprPlaceholderDefault(default: TAT.Expr, value: TAT.Expr, _, _) =>
          inner(default) ++ inner(value)
        case TAT.ExprPlaceholderSep(sep: TAT.Expr, value: TAT.Expr, _, _) =>
          inner(sep) ++ inner(value)
        // Access an array element at [index]
        case TAT.ExprAt(value, index, _, _) =>
          inner(value) ++ inner(index)
        // conditional:
        case TAT.ExprIfThenElse(cond, tBranch, fBranch, _, _) =>
          inner(cond) ++ inner(tBranch) ++ inner(fBranch)
        // Apply a standard library function to arguments.
        case TAT.ExprApply(_, funcWdlType, elements, _, _) =>
          // the function parameters might be optional even if the arguments are not
          def maybeMakeOptional(
              inputs: Vector[WdlInputRef],
              argType: T
          ): Vector[WdlInputRef] = {
            if (TypeUtils.isOptional(argType)) {
              inputs.map(ref => ref.copy(kind = InputKind.Optional))
            } else {
              inputs
            }
          }

          funcWdlType match {
            case _: T_Function0 => Vector.empty
            case T_Function1(_, arg0Type, _) =>
              maybeMakeOptional(inner(elements(0)), arg0Type)
            case T_Function2(_, arg0Type, arg1Type, _) =>
              Vector(
                  maybeMakeOptional(inner(elements(0)), arg0Type),
                  maybeMakeOptional(inner(elements(1)), arg1Type)
              ).flatten
            case T_Function3(_, arg0Type, arg1Type, arg2Type, _) =>
              Vector(
                  maybeMakeOptional(inner(elements(0)), arg0Type),
                  maybeMakeOptional(inner(elements(1)), arg1Type),
                  maybeMakeOptional(inner(elements(2)), arg2Type)
              ).flatten
          }
        // Access a field of a LHS expression that evaluates to a call, pair or object.
        // In the case of a call, we want the fully-qualified name of the field, because
        // the input parameter of the application generated by the compiler will
        // have a fully-qualified name, and that is what we'll need to look up
        // in the evaluation context. In other cases (pair, struct, object), we want a
        // reference to the LHS object, not directly to the field.
        case TAT.ExprGetName(TAT.ExprIdentifier(callId, T_Call(_, output), _),
                             fieldName,
                             wdlType,
                             _) if output.contains(fieldName) =>
          val kind = if (TypeUtils.isOptional(wdlType)) {
            InputKind.Optional
          } else {
            InputKind.Required
          }
          Vector(WdlInputRef(Vector(callId), Some(fieldName), wdlType, kind))
        case TAT.ExprGetName(expr, fieldName, wdlType, _) =>
          // throw an exception if the reference is not valid
          TypeUtils.unwrapOptional(expr.wdlType) match {
            case _: T_Pair if Set("left", "right").contains(fieldName) => ()
            case T_Struct(_, members) if members.contains(fieldName)   => ()
            case T_Object | T_Any                                      => ()
            case _ =>
              throw new Exception(
                  s"Unhandled ExprGetName construction ${TypeUtils.prettyFormatExpr(expr)}"
              )
          }
          // get the LHS inputs, but update optional based on the actual type of this reference
          val kind = if (TypeUtils.isOptional(wdlType)) InputKind.Optional else InputKind.Required
          inner(expr) match {
            case Vector(WdlInputRef(lhsName, Some(field), _, _)) if withField =>
              Vector(WdlInputRef(lhsName :+ field, Some(fieldName), wdlType, kind))
            case Vector(WdlInputRef(lhsName, None, _, _)) if withField =>
              Vector(WdlInputRef(lhsName, Some(fieldName), wdlType, kind))
            case v if withField && v.nonEmpty =>
              throw new Exception(
                  s"cannot add field name because multiple inputs are required to evaluate LHS ${expr}"
              )
            case v =>
              v.map(ref => ref.copy(kind = kind))
          }
        case other =>
          throw new Exception(s"Unhandled expression ${other}")
      }
    }
    inner(expr)
  }

  /**
    * Get all inputs and outputs for a block of statements.
    * @param elements the block elements
    * @return
    */
  def getClosureInputsAndOutputs(
      elements: Vector[TAT.WorkflowElement],
      withField: Boolean
  ): (Map[String, (T, InputKind.InputKind)], Map[String, TAT.OutputParameter]) = {
    def getOutputs(
        innerElements: Vector[TAT.WorkflowElement]
    ): Vector[TAT.OutputParameter] = {
      innerElements.flatMap {
        case TAT.PrivateVariable(name, wdlType, expr, loc) =>
          Vector(TAT.OutputParameter(name, wdlType, expr, loc))
        case call: TAT.Call =>
          call.callee.output.map {
            case (name, wdlType) =>
              val fqn = s"${call.actualName}.${name}"
              TAT.OutputParameter(
                  fqn,
                  wdlType,
                  TAT.ExprIdentifier(fqn, wdlType, call.loc),
                  call.loc
              )
          }.toVector
        case cond: TAT.Conditional =>
          getOutputs(cond.body).map { out =>
            out.copy(wdlType = TypeUtils.ensureOptional(out.wdlType))
          }
        case scatter: TAT.Scatter =>
          // make outputs arrays, remove the collection iteration variable
          val nonEmptyOutputArray = scatter.expr.wdlType match {
            case T_Array(_, nonEmpty) => nonEmpty
            case _ =>
              throw new Exception(
                  s"scatter expression type ${scatter.expr.wdlType} not an array"
              )
          }
          getOutputs(scatter.body).collect {
            case out: TAT.OutputParameter if out.name != scatter.identifier =>
              out.copy(wdlType = T_Array(out.wdlType, nonEmpty = nonEmptyOutputArray))
          }
      }
    }

    def getInputs(
        innerElements: Vector[TAT.WorkflowElement],
        innerWithField: Boolean
    ): Vector[WdlInputRef] = {
      innerElements.flatMap {
        case v: TAT.PrivateVariable =>
          getExpressionInputs(v.expr, innerWithField)
        case call: TAT.Call =>
          call.callee.input.flatMap {
            case (name: String, (_: T, optional: Boolean)) =>
              call.inputs
                .get(name)
                .map { expr =>
                  val exprInputs = getExpressionInputs(expr, withField = true)
                  if (optional) {
                    exprInputs.map { ref =>
                      ref.copy(kind = InputKind.Optional)
                    }
                  } else {
                    exprInputs
                  }
                }
                .getOrElse(Vector.empty)
          }.toVector
        case cond: TAT.Conditional =>
          // get inputs for the conditional expression
          val exprInputs = getExpressionInputs(cond.expr, innerWithField)
          // recurse into body of conditional
          val bodyInputs = getInputs(cond.body, innerWithField = true)
          exprInputs ++ bodyInputs
        case scatter: TAT.Scatter =>
          // get inputs for the scatter expression
          val exprInputs = getExpressionInputs(scatter.expr, innerWithField)
          // recurse into body of the scatter
          // if the scatter variable is referenced, ensure its kind is
          // 'Computed' so it doesn't become a required input
          val scatterIdentifierRegexp = s"${scatter.identifier}([.\\[].+)?".r
          val bodyInputs = getInputs(scatter.body, innerWithField = true).map {
            case ref if scatterIdentifierRegexp.matches(ref.identifier) =>
              ref.copy(kind = InputKind.Computed)
            case ref => ref
          }
          exprInputs ++ bodyInputs
      }
    }

    // first convert outputs - we need to do this prior to
    // inputs because WDL allows forward references, and we
    // need to be able to distinguish block inputs from
    // variables that are defined within the block
    val outputs = getOutputs(elements).groupBy(_.name).map {
      case (name, outputs) if outputs.size == 1 => name -> outputs.head
      case (name, outputs) if Set(outputs).size > 1 =>
        throw new Exception(s"multiple outputs defined with the name ${name}: ${outputs}")
    }

    // now convert the inputs, filter out those that are in outputs,
    // and check for collisions
    val inputs = getInputs(elements, withField)
      .filterNot(i => i.nameIter.exists(outputs.contains))
      .groupBy(_.fullyQualifiedName)
      .map {
        case (fqn, refs) if refs.toSet.size == 1 =>
          fqn -> (refs.head.wdlType, refs.head.kind)
        case (fqn, refs) =>
          // there are multiple references to the same variable from different kinds of
          // input - sort by InputKind, pick the one with the highest priority (lowest
          // value), and make sure if there are multiple with the same kind that they're
          // all of the same type
          val priorityRefs = refs.groupBy(_.kind).toVector.sortWith(_._1 < _._1).head._2
          if (priorityRefs.map(_.wdlType).toSet.size > 1) {
            throw new Exception(
                s"multiple references to the same paramter with different types: ${priorityRefs}"
            )
          }
          fqn -> (priorityRefs.head.wdlType, priorityRefs.head.kind)
      }

    (inputs, outputs)
  }

  /**
    * We are building an applet for the output section of a workflow. The outputs have
    * expressions, and we need to figure out which variables they refer to. This will
    * allow the calculations to proceeed inside a stand alone applet.
    * @param outputs output definitions
    * @return
    */
  def getOutputClosure(outputs: Vector[TAT.OutputParameter]): Map[String, T] = {
    // create inputs from all the expressions that go into outputs
    outputs
      .flatMap {
        case TAT.OutputParameter(_, _, expr, _) => Vector(expr)
      }
      .flatMap(e => getExpressionInputs(e, withField = false))
      .groupBy(_.fullyQualifiedName)
      .map {
        // if there are multiple references to the same parameter, make sure the
        // types are the same
        case (fqn, refs) if refs.toSet.size == 1 => fqn -> refs.head.wdlType
        case (fqn, refs) =>
          val priorityRef = refs.groupBy(_.kind).toVector.sortWith(_._1 < _._1).head._2
          if (priorityRef.map(_.wdlType).toSet.size > 1) {
            throw new Exception(
                s"multiple references to the same paramter with different types: ${priorityRef}"
            )
          }
          fqn -> priorityRef.head.wdlType
      }
  }

  def prettyFormatElement(element: TAT.WorkflowElement, indent: String = "    "): String = {
    element match {
      case TAT.Scatter(varName, expr, body, _) =>
        val collection = TypeUtils.prettyFormatExpr(expr)
        val innerBlock = body
          .map { innerElement =>
            prettyFormatElement(innerElement, indent + "  ")
          }
          .mkString("\n")
        s"""|${indent}scatter (${varName} in ${collection}) {
            |${innerBlock}
            |${indent}}
            |""".stripMargin

      case TAT.Conditional(expr, body, _) =>
        val innerBlock =
          body
            .map { innerElement =>
              prettyFormatElement(innerElement, indent + "  ")
            }
            .mkString("\n")
        s"""|${indent}if (${TypeUtils.prettyFormatExpr(expr)}) {
            |${innerBlock}
            |${indent}}
            |""".stripMargin

      case call: TAT.Call =>
        val inputNames = call.inputs
          .map {
            case (key, expr) =>
              s"${key} = ${TypeUtils.prettyFormatExpr(expr)}"
          }
          .mkString(",")
        val inputs =
          if (inputNames.isEmpty) ""
          else s"{ input: ${inputNames} }"
        call.alias match {
          case None =>
            s"${indent}call ${call.fullyQualifiedName} ${inputs}"
          case Some(al) =>
            s"${indent}call ${call.fullyQualifiedName} as ${al} ${inputs}"
        }

      case TAT.PrivateVariable(name, wdlType, expr, _) =>
        s"${indent} ${TypeUtils.prettyFormatType(wdlType)} ${name} = ${TypeUtils.prettyFormatExpr(expr)}"
    }
  }

  def prettyFormatEnv(env: Map[String, (T, V)], indent: String = "  "): String = {
    env
      .map {
        case (name, (t, v)) =>
          s"${indent}${name}: ${TypeUtils.prettyFormatType(t)} ${EvalUtils.prettyFormat(v)}"
      }
      .mkString("\n")
  }

  def prettyFormatValues(env: Map[String, V], indent: String = "  "): String = {
    env
      .map {
        case (name, v) =>
          s"${indent}${name}: ${EvalUtils.prettyFormat(v)}"
      }
      .mkString("\n")
  }
}

case class IrToWdlValueBindings(
    values: Map[String, Value],
    allowNonstandardCoercions: Boolean = false,
    private var cache: Map[String, V] = Map.empty
) extends VBindings {
  override protected val elementType: String = "value"

  override def contains(name: String): Boolean = values.contains(name)

  override def keySet: Set[String] = values.keySet

  private def resolve(name: String): Unit = {
    if (!cache.contains(name) && values.contains(name)) {
      cache += (name -> WdlUtils.fromIRValue(values(name), Some(name)))
    }
  }

  override def apply(name: String): V = {
    resolve(name)
    cache(name)
  }

  override def get(name: String): Option[V] = {
    resolve(name)
    cache.get(name)
  }

  override lazy val toMap: Map[String, V] = {
    values.keySet.diff(cache.keySet).foreach(resolve)
    cache
  }

  override protected def copyFrom(values: Map[String, V]): IrToWdlValueBindings = {
    copy(cache = cache ++ values)
  }
}

case class ValueMap(values: Map[String, V]) {
  def contains(id: String): Boolean = values.contains(id)

  def get(id: String, wdlTypes: Vector[T] = Vector.empty): Option[V] = {
    (values.get(id), wdlTypes) match {
      case (None, _)         => None
      case (value, Vector()) => value
      case (Some(value), _)  => Some(Coercion.coerceToFirst(wdlTypes, value))
    }
  }
}

object ValueMap {
  lazy val empty: ValueMap = ValueMap(Map.empty)
}
