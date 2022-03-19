package dx.core.languages.wdl

import dx.core.ir.{Application, Callable, ExecutableKindApplet, ExecutableKindNative, SourceCode}
import dx.util.{Logger, StringFileNode}
import spray.json.JsValue
import wdlTools.eval.{DefaultEvalPaths, IoSupport, WdlValues}
import wdlTools.generators.code.{Utils => GeneratorUtils}
import wdlTools.syntax.{CommentMap, Quoting, SourceLocation, WdlVersion}
import wdlTools.types.{GraphUtils, TypeGraph, WdlTypes, TypedAbstractSyntax => TAT}

import scala.collection.immutable.SeqMap

case class WdlDocumentSource(doc: TAT.Document, versionSupport: VersionSupport) extends SourceCode {
  override val language: String = "wdl"

  override def toString: String = versionSupport.generateDocument(doc)

  override def optionsToJson: JsValue = versionSupport.wdlOptions.toJson

  def getDocContents: String = {
    doc.source.readLines.mkString
  }
}

case class WdlWorkflowSource(workflow: TAT.Workflow, versionSupport: VersionSupport)
    extends SourceCode {
  override val language: String = "wdl"

  override def toString: String = versionSupport.generateElement(workflow)

  override def optionsToJson: JsValue = versionSupport.wdlOptions.toJson

  def getDocContents: String = {
    workflow.source.readLines.mkString
  }
}

case class CodeGenerator(typeAliases: Map[String, WdlTypes.T_Struct],
                         wdlVersion: WdlVersion,
                         logger: Logger = Logger.get) {
  // A self contained WDL workflow
  private val outputWdlVersion: WdlVersion = {
    if (wdlVersion == WdlVersion.Draft_2) {
      logger.warning("Upgrading draft-2 input to version 1.0")
      WdlVersion.V1
    } else {
      wdlVersion
    }
  }

  private lazy val structDefs: Vector[TAT.StructDefinition] = {
    val ordered: Vector[(String, WdlTypes.T_Struct)] = if (typeAliases.size <= 1) {
      typeAliases.toVector
    } else {
      // Determine dependency order
      val typeAliasGraph = TypeGraph.buildFromStructTypes(typeAliases)
      GraphUtils
        .toOrderedVector(typeAliasGraph)
        .map { name =>
          typeAliases(name) match {
            case wdlType: WdlTypes.T_Struct =>
              name -> wdlType
            case other =>
              throw new RuntimeException(s"Unexpected type alias ${other}")
          }
        }
    }
    ordered.map {
      case (name, wdlType) =>
        TAT.StructDefinition(name, wdlType, wdlType.members)(SourceLocation.empty)
    }
  }

  private[wdl] def wdlValueToExpr(value: WdlValues.V): TAT.Expr = {
    def seqToType(vec: Iterable[TAT.Expr]): WdlTypes.T = {
      vec.headOption.map(_.wdlType).getOrElse(WdlTypes.T_Any)
    }

    value match {
      case WdlValues.V_Null => TAT.ValueNull(WdlTypes.T_Any)(SourceLocation.empty)
      case WdlValues.V_Boolean(value) =>
        TAT.ValueBoolean(value, WdlTypes.T_Boolean)(SourceLocation.empty)
      case WdlValues.V_Int(value)   => TAT.ValueInt(value, WdlTypes.T_Int)(SourceLocation.empty)
      case WdlValues.V_Float(value) => TAT.ValueFloat(value, WdlTypes.T_Float)(SourceLocation.empty)
      case WdlValues.V_String(value) if value.contains('"') && value.contains("'") =>
        TAT.ValueString(GeneratorUtils.escape(value), WdlTypes.T_String, quoting = Quoting.Double)(
            SourceLocation.empty
        )
      case WdlValues.V_String(value) if value.contains('"') =>
        TAT.ValueString(value, WdlTypes.T_String, quoting = Quoting.Single)(SourceLocation.empty)
      case WdlValues.V_String(value) =>
        TAT.ValueString(value, WdlTypes.T_String, quoting = Quoting.Double)(SourceLocation.empty)
      case WdlValues.V_File(value) => TAT.ValueFile(value, WdlTypes.T_File)(SourceLocation.empty)
      case WdlValues.V_Directory(value) =>
        TAT.ValueDirectory(value, WdlTypes.T_Directory)(SourceLocation.empty)
      // compound values
      case WdlValues.V_Pair(l, r) =>
        val lExpr = wdlValueToExpr(l)
        val rExpr = wdlValueToExpr(r)
        TAT.ExprPair(lExpr, rExpr, WdlTypes.T_Pair(lExpr.wdlType, rExpr.wdlType))(
            SourceLocation.empty
        )
      case WdlValues.V_Array(value) =>
        val valueExprs = value.map(wdlValueToExpr)
        TAT.ExprArray(valueExprs, seqToType(valueExprs))(SourceLocation.empty)
      case WdlValues.V_Map(value) =>
        val keyExprs = value.keys.map(wdlValueToExpr)
        val valueExprs = value.values.map(wdlValueToExpr)
        TAT.ExprMap(
            keyExprs.zip(valueExprs).to(SeqMap),
            WdlTypes.T_Map(seqToType(keyExprs), seqToType(valueExprs))
        )(SourceLocation.empty)
      case WdlValues.V_Optional(value) => wdlValueToExpr(value)
      case WdlValues.V_Struct(name, members) =>
        val memberExprs: Map[TAT.Expr, TAT.Expr] = members.map {
          case (name, value) =>
            TAT.ValueString(name, WdlTypes.T_String)(SourceLocation.empty) -> wdlValueToExpr(value)
          case other => throw new RuntimeException(s"Unexpected member ${other}")
        }
        val memberTypes = memberExprs.map {
          case (name: TAT.ValueString, value) => name.value -> value.wdlType
          case other                          => throw new RuntimeException(s"Unexpected member ${other}")
        }
        TAT.ExprMap(memberExprs.to(SeqMap), WdlTypes.T_Struct(name, memberTypes.to(SeqMap)))(
            SourceLocation.empty
        )
      case WdlValues.V_Object(members) =>
        val memberExprs = members.map {
          case (name, value) =>
            val key: TAT.Expr = TAT.ValueString(name, WdlTypes.T_String)(SourceLocation.empty)
            key -> wdlValueToExpr(value)
        }
        TAT.ExprObject(memberExprs.to(SeqMap), WdlTypes.T_Object)(SourceLocation.empty)
      case other => throw new Exception(s"Unhandled value ${other}")
    }
  }

  /*
  Create a header for a task/workflow. This is an empty task that includes the input and output
  definitions. It is used to
  (1) allow linking to native DNAx applets (and workflows in the future).
  (2) make a WDL file stand-alone, without imports

  For example, the stub for the Add task:

  task Add {
    input {
      Int a
      Int b
    }
    command <<<
    python -c "print(${a} + ${b})"
    >>>
    output {
      Int result = read_int(stdout())
    }
  }

  is:

  task Add {
    input {
      Int a
      Int b
    }
    command {}
    output {
      Int result
    }
  }
   */
  private def createTaskStub(callable: Callable,
                             native: Option[ExecutableKindNative] = None): TAT.Task = {
    // Sort the inputs by name, so the result will be deterministic.
    val inputs: Vector[TAT.InputParameter] =
      callable.inputVars
        .sortWith(_.name < _.name)
        .map { parameter =>
          val wdlType = WdlUtils.fromIRType(parameter.dxType, typeAliases)
          parameter.defaultValue match {
            case None =>
              TAT.RequiredInputParameter(parameter.name.decoded, wdlType)(SourceLocation.empty)
            case Some(value) =>
              TAT.OverridableInputParameterWithDefault(
                  parameter.name.decoded,
                  wdlType,
                  WdlUtils.irValueToExpr(value)
              )(SourceLocation.empty)
          }
        }
    val outputs: Vector[TAT.OutputParameter] =
      callable.outputVars
        .sortWith(_.name < _.name)
        .map { parameter =>
          val wdlType = WdlUtils.fromIRType(parameter.dxType, typeAliases)
          val defaultVal = WdlUtils.getDefaultValueOfType(wdlType)
          TAT.OutputParameter(parameter.name.decoded, wdlType, defaultVal)(SourceLocation.empty)
        }
    val meta = native.map { n =>
      TAT.MetaSection(
          Vector(
              Some("type" -> TAT.MetaValueString("native", Quoting.Double)(SourceLocation.empty)),
              n.id.map(id => "id" -> TAT.MetaValueString(id, Quoting.Double)(SourceLocation.empty))
          ).flatten.to(SeqMap)
      )(
          SourceLocation.empty
      )
    }
    TAT.Task(
        callable.name,
        WdlTypes.T_Task(
            callable.name,
            inputs
              .map {
                case TAT.RequiredInputParameter(name, wdlType) =>
                  name -> (wdlType, false)
                case other: TAT.InputParameter =>
                  other.name -> (other.wdlType, true)
              }
              .to(SeqMap),
            outputs.map(d => d.name -> d.wdlType).to(SeqMap),
            None
        ),
        inputs,
        outputs,
        TAT.CommandSection(Vector.empty)(SourceLocation.empty),
        Vector.empty,
        meta,
        None,
        None,
        None
    )(
        SourceLocation.empty
    )
  }

  /**
    * Generate a WDL stub for a DNAnexus applet.
    * @param id the applet ID
    * @param appletName the applet name
    * @param inputSpec the applet inputs
    * @param outputSpec the applet outputs
    * @return an AST.Task
    */
  def createAppletStub(id: String,
                       appletName: String,
                       inputSpec: Map[String, WdlTypes.T],
                       outputSpec: Map[String, WdlTypes.T]): TAT.Task = {
    TAT.Task(
        appletName,
        WdlTypes.T_Task(appletName,
                        inputSpec
                          .map {
                            case (name, wdlType) => name -> (wdlType, false)
                          }
                          .to(SeqMap),
                        outputSpec.to(SeqMap),
                        None),
        inputSpec.map {
          case (name, wdlType) => TAT.RequiredInputParameter(name, wdlType)(SourceLocation.empty)
        }.toVector,
        outputSpec.map {
          case (name, wdlType) =>
            val expr = WdlUtils.getDefaultValueOfType(wdlType)
            TAT.OutputParameter(name, wdlType, expr)(SourceLocation.empty)
        }.toVector,
        TAT.CommandSection(Vector.empty)(SourceLocation.empty),
        Vector.empty,
        Some(
            TAT.MetaSection(
                SeqMap(
                    "type" -> TAT.MetaValueString("native", Quoting.Double)(SourceLocation.empty),
                    "id" -> TAT.MetaValueString(id, Quoting.Double)(SourceLocation.empty)
                )
            )(
                SourceLocation.empty
            )
        ),
        parameterMeta = None,
        runtime = None,
        hints = None
    )(
        SourceLocation.empty
    )
  }

  def createStandAloneTask(task: TAT.Task): TAT.Document = {
    val ioSupp = IoSupport(DefaultEvalPaths.empty)
    val srcString = ioSupp.readFilePosition(task.loc.source.toString, task.loc)
    TAT.Document(StringFileNode(contents = srcString),
                 TAT.Version(outputWdlVersion)(SourceLocation.empty),
                 structDefs :+ task,
                 None,
                 CommentMap.empty)(
        SourceLocation.empty
    )
  }

  /*
   A workflow can import other libraries:

   import "library.wdl" as lib
   workflow foo {
     call lib.Multiply as mul { ... }
     call lib.Add { ... }
     call lib.Nice as nice { ... }
     call lib.Hello
   }

   rewrite the workflow, and remove the calls to external libraries.

   workflow foo {
     call Multiply as mul { ... }
     call Add { ... }
     call Nice as nice { ... }
     call Hello
   }
   */
  private def unqualifyCallNames(body: Vector[TAT.WorkflowElement]): Vector[TAT.WorkflowElement] = {
    body.map {
      case call: TAT.Call =>
        call.copy(fullyQualifiedName = call.unqualifiedName)(call.loc)
      case scatter: TAT.Scatter =>
        scatter.copy(body = unqualifyCallNames(scatter.body))(scatter.loc)
      case cond: TAT.Conditional =>
        cond.copy(body = unqualifyCallNames(cond.body))(cond.loc)
      case other => other
    }
  }

  /**
    * A workflow must have definitions for all the tasks it calls. However, a scatter calls tasks
    * that are missing from the WDL file we generate. To ameliorate this, we add stubs for called
    * tasks. The generated tasks are named by their unqualified names, not their fully-qualified
    * names. This works because the WDL workflow must be "flattenable".
    * @param wf the workflow
    * @param callables the callables to add to the workflow
    * @return
    */
  def standAloneWorkflow(wf: TAT.Workflow, callables: Vector[Callable]): TAT.Document = {
    val tasks: Vector[TAT.Task] =
      callables
        .foldLeft(Map.empty[String, TAT.Task]) {
          case (accu, callable) =>
            if (accu contains callable.name) {
              // we have already created a stub for this call
              accu
            } else {
              val stub: TAT.Task = callable match {
                case Application(_,
                                 _,
                                 _,
                                 _,
                                 _,
                                 ExecutableKindApplet,
                                 WdlDocumentSource(doc, _),
                                 _,
                                 _,
                                 _,
                                 _,
                                 _) =>
                  // This is a task, include its source instead of a header.
                  val tasks = doc.elements.collect {
                    case t: TAT.Task => t
                  }
                  assert(tasks.size == 1)
                  tasks.head
                case Application(_, _, _, _, _, native: ExecutableKindNative, _, _, _, _, _, _) =>
                  // no existing stub, create it - specify whether the target app(let) is native
                  createTaskStub(callable, native = Some(native))
                case _ =>
                  // no existing stub, create it
                  createTaskStub(callable)
              }
              accu + (callable.name -> stub)
            }
        }
        // sort the task order by name, so the generated code will be deterministic
        .toVector
        .sortWith(_._1 < _._1)
        .map { case (_, task) => task }

    val wfWithoutImportCalls = wf.copy(body = unqualifyCallNames(wf.body))(wf.loc)
    val ioSupp = IoSupport(DefaultEvalPaths.empty)
    val srcString = ioSupp.readFilePosition(wf.loc.source.toString, wf.loc)
    TAT.Document(
        StringFileNode(contents = srcString),
        TAT.Version(outputWdlVersion)(SourceLocation.empty),
        structDefs ++ tasks,
        Some(wfWithoutImportCalls),
        CommentMap.empty
    )(SourceLocation.empty)
  }
}
