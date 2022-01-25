package dx.translator

import java.nio.file.{Path, Paths}
import dx.Tags.EdgeTest
import dx.api._
import dx.core.Constants
import dx.core.CliUtils.{Failure, UnsuccessfulTermination}
import dx.core.ir.{Parameter, _}
import dx.core.ir.RunSpec._
import dx.core.ir.Type._
import dx.core.ir.Value._
import dx.core.languages.cwl.{CwlDxName, TargetParam}
import dx.core.languages.wdl.{WdlDocumentSource, WdlDxName}
import dx.translator.CallableAttributes._
import dx.translator.ParameterAttributes._
import dx.util.Logger
import dxCompiler.Main
import dxCompiler.Main.SuccessfulCompileIR
import org.scalatest.Inside._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdlTools.generators.code.WdlGenerator

import scala.collection.immutable.SeqMap

// These tests involve compilation -without- access to the platform.
//
class TranslatorTest extends AnyFlatSpec with Matchers {
  private val dxApi = DxApi()(Logger.Quiet)

  private def pathFromBasename(dir: String, basename: String): Path = {
    Paths.get(getClass.getResource(s"/${dir}/${basename}").getPath)
  }

  private val dxProjectId = dxApi.currentProjectId.get

  // task compilation
  private val cFlags =
    List("--compileMode", "ir", "-quiet", "--locked", "--project", dxProjectId)

  private val cFlagsUnlocked =
    List("--compileMode", "ir", "-quiet", "--project", dxProjectId)

  val dbgFlags = List("--compileMode",
                      "ir",
                      "--verbose",
                      "--verboseKey",
                      "GenerateIR",
                      "--locked",
                      "--project",
                      dxProjectId)

  private def getApplicationByName(name: String, bundle: Bundle): Application =
    bundle.allCallables(name) match {
      case a: Application => a
      case _              => throw new Exception(s"${name} is not an applet")
    }

  it should "IR compile a single WDL task" in {
    val path = pathFromBasename("compiler", "add.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a single WDL task with dynamic instance type selection" in {
    def compile(instanceTypeSelection: String): InstanceType = {
      val path = pathFromBasename("compiler", "add.wdl")
      val args = path.toString :: "-instanceTypeSelection" :: instanceTypeSelection :: cFlags
      val bundle = Main.compile(args.toVector) match {
        case SuccessfulCompileIR(bundle) => bundle
        case other =>
          throw new Exception(s"expected succss not ${other}")
      }
      val applet = bundle.primaryCallable match {
        case Some(applet: Application) => applet
        case other =>
          throw new Exception(s"expected primary callable to be an applet not ${other}")
      }
      applet.instanceType
    }
    compile("static") shouldBe StaticInstanceType(
        InstanceTypeRequest(None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            Some(ExecutionEnvironment("Ubuntu", "20.04", Vector("0"))))
    )
    compile("dynamic") shouldBe DynamicInstanceType
  }

  it should "IR compile a task with docker" in {
    val path = pathFromBasename("compiler", "BroadGenomicsDocker.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  // workflow compilation
  it should "IR compile a linear WDL workflow without expressions" in {
    val path = pathFromBasename("compiler", "wf_linear_no_expr.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a linear WDL workflow" in {
    val path = pathFromBasename("compiler", "wf_linear.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile unlocked workflow" in {
    val path = pathFromBasename("compiler", "wf_linear.wdl")
    val args = path.toString :: cFlagsUnlocked
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a non trivial linear workflow with variable coercions" in {
    val path = pathFromBasename("compiler", "cast.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a workflow with two consecutive calls" in {
    val path = pathFromBasename("compiler", "strings.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a workflow with a scatter without a call" in {
    val path = pathFromBasename("compiler", "scatter_no_call.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile optionals" in {
    val path = pathFromBasename("compiler", "optionals.wdl")
    val args = path.toString :: cFlags
    //                :: "--verbose"
    //                :: "--verboseKey" :: "GenerateIR"
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "support imports" in {
    val path = pathFromBasename("compiler", "check_imports.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "IR compile a draft2 workflow" in {
    val path = pathFromBasename("draft2", "shapes.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "expressions in an output block" in {
    val path = pathFromBasename("compiler", "expr_output_block.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  /*  ignore should "scatters over maps" in {
    val path = pathFromBasename("compiler", "dict2.wdl")
    val args =         path.toString :: cFlags
Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }*/

  it should "skip missing optional arguments" in {
    val path = pathFromBasename("util", "missing_inputs_to_direct_call.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "handle calling subworkflows" in {
    val path = pathFromBasename("subworkflows", "trains.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(irwf) => irwf
      case _                         => throw new Exception("unexpected")
    }
    val primaryWf: Workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case _                  => throw new Exception("unexpected")
    }
    primaryWf.stages.size shouldBe 2
  }

  it should "compile a sub-block with several calls" in {
    val path = pathFromBasename("compiler", "subblock_several_calls.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "missing workflow inputs" in {
    val path = pathFromBasename("input_file", "missing_args.wdl")
    val args = path.toString :: List("--compileMode", "ir", "--quiet", "--project", dxProjectId)
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  // Nested blocks
  it should "compile two level nested workflow" in {
    val path = pathFromBasename("nested", "two_levels.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "handle passing closure arguments to nested blocks" in {
    val path = pathFromBasename("nested", "param_passing.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "compile a workflow calling a subworkflow as a direct call" in {
    val path = pathFromBasename("draft2", "movies.wdl")
    val args = path.toString :: cFlags
    val bundle: Bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        Logger.error(other.toString)
        throw new Exception(s"Failed to compile ${path}")
    }
    val wf: Workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) =>
        wf
      case _ => throw new Exception("bad value in bundle")
    }
    val stage = wf.stages.head
    stage.description shouldBe "review"
  }

  it should "compile a workflow calling a subworkflow as a direct call with 2.0 version" in {
    val path = pathFromBasename("v2", "movies.wdl")
    val args = path.toString :: cFlags
    val bundle: Bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        Logger.error(other.toString)
        throw new Exception(s"Failed to compile ${path}")
    }
    val wf: Workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) =>
        wf
      case _ => throw new Exception("bad value in bundle")
    }
    val stage = wf.stages.head
    stage.description shouldBe "review"
  }

  it should "compile a workflow calling a subworkflow with native DNANexus applet as a direct call with 2.0 version" in {
    val path = pathFromBasename("v2", "call_dnanexus_applet.wdl")
    val args = path.toString :: cFlags
    val bundle: Bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        Logger.error(other.toString)
        throw new Exception(s"Failed to compile ${path}")
    }
    val wf: Workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) =>
        wf
      case _ => throw new Exception("bad value in bundle")
    }
    wf.stages.size shouldBe 2
    wf.stages(0).description shouldBe "native_sum_012"
    wf.stages(1).description shouldBe "native_sum_wf"
  }

  it should "three nesting levels" in {
    val path = pathFromBasename("nested", "three_levels.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val primary: Callable = bundle.primaryCallable.get
    val wf = primary match {
      case wf: Workflow => wf
      case _            => throw new Exception("unexpected")
    }

    wf.stages.size shouldBe 1

    val level2 = bundle.allCallables(wf.name)
    level2 shouldBe a[Workflow]
    val wfLevel2 = level2.asInstanceOf[Workflow]
    wfLevel2.stages.size shouldBe 1
  }

  it should "four nesting levels" in {
    val path = pathFromBasename("nested", "four_levels.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a workflow with nested scatter" in {
    val path = pathFromBasename("nested", "nested_scatter.wdl")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other                       => throw new Exception(s"unexpected compile result ${other}")
    }
    val outerScatter = bundle.allCallables("nested_scatter_frag_stage-3")
    outerScatter.inputVars should contain theSameElementsAs Vector(
        Parameter(WdlDxName.fromSourceName("x"), TInt),
        Parameter(WdlDxName.fromSourceName("ints1"), TArray(TInt)),
        Parameter(WdlDxName.fromSourceName("ints2"), TArray(TInt))
    )
    val innerScatter = bundle.allCallables("nested_scatter_block_0_0")
    innerScatter.inputVars should contain theSameElementsAs Vector(
        Parameter(WdlDxName.fromSourceName("x"), TInt),
        Parameter(WdlDxName.fromSourceName("y"), TInt),
        Parameter(WdlDxName.fromSourceName("ints1"), TArray(TInt))
    )
  }

  // Check parameter_meta `pattern` keyword
  it should "recognize pattern in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "pattern_params.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("pattern_params_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            Type.TFile,
            None,
            Vector(
                HelpAttribute("The input file to be searched"),
                PatternsAttribute(PatternsArray(Vector("*.txt", "*.tsv"))),
                GroupAttribute("Common"),
                LabelAttribute("Input file")
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                HelpAttribute("The pattern to use to search in_file"),
                GroupAttribute("Common"),
                LabelAttribute("Search pattern")
            )
        )
    )
    cgrepApplication.outputs.iterator sameElements Vector(
        Parameter(WdlDxName.fromSourceName("count"), TInt, None, Vector.empty),
        Parameter(
            WdlDxName.fromSourceName("out_file"),
            TFile,
            None,
            Vector(
                PatternsAttribute(PatternsArray(Vector("*.txt", "*.tsv"))),
                GroupAttribute("Common"),
                LabelAttribute("Output file")
            )
        )
    )
  }

  // Check parameter_meta `pattern` keyword
  it should "recognize pattern object in parameters_obj_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "pattern_obj_params.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("pattern_params_obj_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                HelpAttribute("The input file to be searched"),
                PatternsAttribute(
                    PatternsObject(
                        Vector("*.txt", "*.tsv"),
                        Some("file"),
                        Vector("foo", "bar")
                    )
                ),
                GroupAttribute("Common"),
                LabelAttribute("Input file")
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                HelpAttribute("The pattern to use to search in_file"),
                GroupAttribute("Common"),
                LabelAttribute("Search pattern")
            )
        )
    )
    cgrepApplication.outputs.iterator sameElements Vector(
        Parameter(WdlDxName.fromSourceName("count"), TInt, None, Vector.empty),
        Parameter(
            WdlDxName.fromSourceName("out_file"),
            TFile,
            None,
            Vector(
                PatternsAttribute(PatternsArray(Vector("*.txt", "*.tsv"))),
                GroupAttribute("Common"),
                LabelAttribute("Input file")
            )
        )
    )
  }

  // Check parameter_meta `choices` keyword
  it should "recognize choices in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "choice_values.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("choice_values_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                ChoicesAttribute(
                    Vector(
                        FileChoice(
                            name = None,
                            value = "dx://file-Fg5PgBQ0ffP7B8bg3xqB115G"
                        ),
                        FileChoice(
                            name = None,
                            value = "dx://file-Fg5PgBj0ffPP0Jjv3zfv0yxq"
                        )
                    )
                )
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                ChoicesAttribute(
                    Vector(
                        SimpleChoice(value = VString("A")),
                        SimpleChoice(value = VString("B"))
                    )
                )
            )
        )
    )
  }

  // Check parameter_meta `choices` keyword with annotated values
  it should "recognize annotated choices in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "choice_obj_values.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("choice_values_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                ChoicesAttribute(
                    Vector(
                        FileChoice(
                            name = Some("file1"),
                            value = "dx://file-Fg5PgBQ0ffP7B8bg3xqB115G"
                        ),
                        FileChoice(
                            name = Some("file2"),
                            value = "dx://file-Fg5PgBj0ffPP0Jjv3zfv0yxq"
                        )
                    )
                )
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                ChoicesAttribute(
                    Vector(
                        SimpleChoice(value = VString("A")),
                        SimpleChoice(value = VString("B"))
                    )
                )
            )
        )
    )
  }

  // Check parameter_meta `suggestion` keyword fails when there is a type mismatch
  it should "throw exception when choice types don't match parameter types" in {
    val path = pathFromBasename("compiler", "choices_type_mismatch.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[UnsuccessfulTermination]
    // TODO: make assertion about exception message
  }

  // Check parameter_meta `suggestions` keyword
  it should "recognize suggestions in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "suggestion_values.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("suggestion_values_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                SuggestionsAttribute(
                    Vector(
                        FileSuggestion(
                            name = None,
                            value = Some("dx://file-Fg5PgBQ0ffP7B8bg3xqB115G"),
                            project = None,
                            path = None
                        ),
                        FileSuggestion(
                            name = None,
                            value = Some("dx://file-Fg5PgBj0ffPP0Jjv3zfv0yxq"),
                            project = None,
                            path = None
                        )
                    )
                )
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                SuggestionsAttribute(
                    Vector(
                        SimpleSuggestion(value = VString("A")),
                        SimpleSuggestion(value = VString("B"))
                    )
                )
            )
        )
    )
  }

  // Check parameter_meta `suggestions` keyword with annotated values
  it should "recognize annotated suggestions in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "suggestion_obj_values.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("suggestion_values_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                SuggestionsAttribute(
                    Vector(
                        FileSuggestion(
                            name = Some("file1"),
                            value = Some("dx://file-Fg5PgBQ0ffP7B8bg3xqB115G"),
                            project = None,
                            path = None
                        ),
                        FileSuggestion(
                            name = Some("file2"),
                            value = None,
                            project = Some("project-FGpfqjQ0ffPF1Q106JYP2j3v"),
                            path = Some("/test_data/f2.txt.gz")
                        )
                    )
                )
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                SuggestionsAttribute(
                    Vector(
                        SimpleSuggestion(value = VString("A")),
                        SimpleSuggestion(value = VString("B"))
                    )
                )
            )
        )
    )
  }

  // Check parameter_meta `suggestions` keyword fails when there is a parameter mismatch
  it should "throw exception when suggestion types don't match parameter types" in {
    val path = pathFromBasename("compiler", "suggestions_type_mismatch.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[UnsuccessfulTermination]
    // TODO: make assertion about exception message
  }

  // Check parameter_meta `suggestions` keyword fails when there is a missing keyword
  it should "throw exception when file suggestion is missing a keyword" in {
    val path = pathFromBasename("compiler", "suggestions_missing_arg.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[UnsuccessfulTermination]
    // TODO: make assertion about exception message
  }

  // Check parameter_meta `dx_type` keyword
  it should "recognize dx_type in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "add_dx_type.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("add_dx_type", bundle)
    cgrepApplication.inputs shouldBe Vector(
        Parameter(
            WdlDxName.fromSourceName("a"),
            TFile,
            None,
            Vector(
                TypeAttribute(DxConstraintString("fastq"))
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("b"),
            TFile,
            None,
            Vector(
                TypeAttribute(
                    DxConstraintBool(
                        DxConstraintOper.And,
                        DxConstraintArray(
                            DxConstraintString("fastq"),
                            DxConstraintBool(
                                DxConstraintOper.Or,
                                DxConstraintArray(
                                    DxConstraintString("Read1"),
                                    DxConstraintString("Read2")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
  }

  // Check parameter_meta `dx_type` keyword fails when specified for a non-file parameter
  it should "throw exception when dx_type is used on non-file parameter" in {
    val path = pathFromBasename("compiler", "dx_type_nonfile.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[UnsuccessfulTermination]
    // TODO: make assertion about exception message
  }

  // Check parameter_meta `default` keyword
  it should "recognize default in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "add_default.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("add_default", bundle)
    cgrepApplication.inputs shouldBe Vector(
        Parameter(
            WdlDxName.fromSourceName("a"),
            TInt,
            Some(VInt(1)),
            Vector.empty
        ),
        Parameter(
            WdlDxName.fromSourceName("b"),
            TOptional(TInt),
            None,
            Vector(DefaultAttribute(VInt(2)))
        )
    )
  }

  // Check parameter_meta `default` keyword fails when there is a type mismatch
  it should "throw exception when default types don't match parameter types" in {
    val path = pathFromBasename("compiler", "default_type_mismatch.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[UnsuccessfulTermination]
    // TODO: make assertion about exception message
  }

  it should "recognize help in parameters_meta via Parameter for input Parameters" in {
    val path = pathFromBasename("compiler", "help_input_params.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("help_input_params_cgrep", bundle)
    cgrepApplication.inputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("s"),
            TString,
            None,
            Vector(
                HelpAttribute("This is help for s")
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("in_file"),
            TFile,
            None,
            Vector(
                HelpAttribute("The input file to be searched"),
                GroupAttribute("Common"),
                LabelAttribute("Input file")
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("pattern"),
            TString,
            None,
            Vector(
                HelpAttribute("The pattern to use to search in_file"),
                GroupAttribute("Common"),
                LabelAttribute("Search pattern")
            )
        )
    )
  }

  // This is actually more of a test to confirm that symbols that are not input
  // variables are ignored. WDL doesn't include a paramMeta member for the output
  // var class anyways, so it's basically impossible for this to happen
  it should "ignore help in parameters_meta via Parameter for output Parameters" in {
    val path = pathFromBasename("compiler", "help_output_params.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("help_output_params_cgrep", bundle)
    cgrepApplication.outputs.iterator sameElements Vector(
        Parameter(
            WdlDxName.fromSourceName("count"),
            TInt,
            None,
            Vector.empty
        )
    )
  }

  it should "recognize app metadata" in {
    val path = pathFromBasename("compiler", "add_app_meta.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("add", bundle)
    cgrepApplication.attributes.iterator sameElements
      Vector(
          DeveloperNotesAttribute("Check out my sick bash expression! Three dolla signs!!!"),
          DescriptionAttribute(
              "Adds two int together. This app adds together two integers and returns the sum"
          ),
          TagsAttribute(Vector("add", "ints")),
          OpenSourceAttribute(true),
          VersionAttribute("1.0"),
          PropertiesAttribute(Map("foo" -> "bar")),
          CategoriesAttribute(Vector("Assembly")),
          DetailsAttribute(
              Map(
                  "contactEmail" -> VString("joe@dev.com"),
                  "upstreamVersion" -> VString("1.0"),
                  "upstreamAuthor" -> VString("Joe Developer"),
                  "upstreamUrl" -> VString("https://dev.com/joe"),
                  "upstreamLicenses" -> VArray(
                      Vector(
                          VString("MIT")
                      )
                  ),
                  "whatsNew" -> VArray(
                      Vector(
                          VHash(
                              SeqMap(
                                  "version" -> VString("1.1"),
                                  "changes" -> VArray(
                                      Vector(
                                          VString("Added parameter --foo"),
                                          VString("Added cowsay easter-egg")
                                      )
                                  )
                              )
                          ),
                          VHash(
                              SeqMap(
                                  "version" -> VString("1.0"),
                                  "changes" -> VArray(
                                      Vector(
                                          VString("Initial version")
                                      )
                                  )
                              )
                          )
                      )
                  )
              )
          ),
          TitleAttribute("Add Ints"),
          TypesAttribute(Vector("Adder"))
      )
  }

  it should "recognize runtime hints" in {
    val path = pathFromBasename("compiler", "add_runtime_hints.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepApplication = getApplicationByName("add_runtime_hints", bundle)
    cgrepApplication.requirements.iterator sameElements
      Vector(
          IgnoreReuseRequirement(true),
          RestartRequirement(
              max = Some(5),
              default = Some(1),
              errors = Map("UnresponsiveWorker" -> 2, "ExecutionError" -> 2)
          ),
          TimeoutRequirement(hours = Some(12), minutes = Some(30)),
          AccessRequirement(network = Vector("*"), developer = Some(true))
      )
  }

  it should "ignore dx_instance_type when evaluating runtime hints" in {
    val path = pathFromBasename("compiler", "instance_type_test.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
  }

  it should "handle an empty workflow" in {
    val path = pathFromBasename("util", "empty_workflow.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle structs" in {
    val path = pathFromBasename("struct", "Person.wdl")
    val args = path.toString :: cFlags
    val retval =
      Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle access to struct member" in {
    val path = pathFromBasename("struct", "struct_deref.wdl")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(ir) => ir
      case other =>
        throw new Exception(s"unexpected result ${other}")
    }
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other =>
        throw new Exception(s"unexpected primary callable ${other}")
    }
    wf.stages.size shouldBe 2
    wf.stages.head.inputs shouldBe Vector(
        StageInputWorkflowLink(
            Parameter(WdlDxName.fromSourceName("sampleStruct"),
                      TSchema("SampleStruct", SeqMap("sample_name" -> TString, "id" -> TInt)),
                      None,
                      Vector())
        )
    )
  }

  it should "recognize that an argument with a default can be omitted at the call site" in {
    val path = pathFromBasename("compiler", "call_level2.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "check for reserved symbols" in {
    val path = pathFromBasename("compiler", "reserved.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    inside(retval) {
      case Failure(_, Some(e)) =>
        e.getMessage should include("using the substring '___'")
    }
  }

  it should "do nested scatters" in {
    val path = pathFromBasename("compiler", "nested_scatter.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle struct imported several times" in {
    val path = pathFromBasename("struct/struct_imported_twice", "file3.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle file constants in a workflow" in {
    val path = pathFromBasename("compiler", "wf_constants.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "respect import flag" in {
    val path = pathFromBasename("compiler/imports", "A.wdl")
    val libraryPath = path.getParent.resolve("lib")
    val args = path.toString :: "--imports" :: libraryPath.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "respect import -p flag" in {
    val path = pathFromBasename("compiler/imports", "A.wdl")
    val libraryPath = path.getParent.resolve("lib")
    val args = path.toString :: "-p" :: libraryPath.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "pass environment between deep stages" in {
    val path = pathFromBasename("compiler", "environment_passing_deep_nesting.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "handle multiple struct definitions" in {
    val path = pathFromBasename("struct/DEVEX-1196-struct-resolution-wrong-order", "file3.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "retain all characters in a WDL task" in {
    val path = pathFromBasename("bugs", "missing_chars_in_task.wdl")
    val args = path.toString :: cFlags
    //                                      :: "--verbose"
    //                                      :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val commandSection =
      """|  command <<<
         |  echo 1 hello world | sed 's/world/wdl/'
         |  echo 2 hello \
         |  world \
         |  | sed 's/world/wdl/'
         |  echo 3 hello \
         |  world | \
         |  sed 's/world/wdl/'
         |  >>>
         |""".stripMargin

    inside(retval) {
      case SuccessfulCompileIR(bundle) =>
        bundle.allCallables.size shouldBe 1
        val (_, callable) = bundle.allCallables.head
        callable shouldBe a[Application]
        val task = callable.asInstanceOf[Application]
        val generator = WdlGenerator()
        val wdlDoc = task.document match {
          case WdlDocumentSource(doc, _) => doc
          case _                         => throw new Exception("expected a WDL document")
        }
        val taskSource = generator.generateDocument(wdlDoc).mkString("\n")
        taskSource should include(commandSection)
    }
  }

  it should "correctly flatten a workflow with imports" in {
    val path = pathFromBasename("compiler", "wf_to_flatten.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "detect a request for GPU" in {
    val path = pathFromBasename("compiler", "GPU.wdl")
    val args = path.toString :: cFlags
    //                                      :: "--verbose"
    //                                      :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    inside(retval) {
      case SuccessfulCompileIR(bundle) =>
        bundle.allCallables.size shouldBe 1
        val (_, callable) = bundle.allCallables.head
        callable shouldBe a[Application]
        val task = callable.asInstanceOf[Application]
        task.instanceType shouldBe StaticInstanceType(
            InstanceTypeRequest(Some("mem3_ssd1_gpu_x8"),
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None,
                                None)
        )
    }
  }

  it should "compile a scatter with a sub-workflow that has an optional argument" in {
    val path = pathFromBasename("compiler", "scatter_subworkflow_with_optional.wdl")
    val args = path.toString :: cFlags
    //                                      :: "--verbose"
    //                                      :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val bundle = retval match {
      case SuccessfulCompileIR(bundle) => bundle
      case _                           => throw new Exception("unexpected")
    }

    val wfs: Vector[Workflow] = bundle.allCallables.flatMap {
      case (_, wf: Workflow) if wf.locked && wf.level == Level.Sub => Some(wf)
      case (_, _)                                                  => None
    }.toVector
    wfs.length shouldBe 1
    val wf = wfs.head

    val samtools = wf.inputs.find {
      case (cVar, _) => cVar.name.decoded == "samtools_memory"
    }
    inside(samtools) {
      /*case Some((cVar, _)) =>
       cVar.wdlType shouldBe (TOptional(TString))*/
      case None => ()
    }
  }
  it should "compile a workflow taking arguments from a Pair" in {
    val path = pathFromBasename("draft2", "pair.wdl")
    val args = path.toString :: cFlags
    //                                      :: "--verbose"
    //                                      :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "pass as subworkflows do not have expression statement in output block" in {
    val path = pathFromBasename("subworkflows", basename = "trains.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  // this is currently failing.
  it should "pass with subworkflows having expression" in {
    val path = pathFromBasename("subworkflows", basename = "ensure_trains.wdl")

    /* ensure_trains workflow
     * trains        workflow
     * check_route   workflow
     * concat        task
     */
    val args = path.toString :: cFlags
    //          :: "--verbose"
    //          :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "recognize workflow metadata" in {
    val path = pathFromBasename("compiler", "wf_meta.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case _                  => throw new Exception("primaryCallable is not a workflow")
    }
    workflow.attributes.iterator sameElements
      Vector(
          DescriptionAttribute("This is a workflow that defines some metadata"),
          TagsAttribute(Vector("foo", "bar")),
          VersionAttribute("1.0"),
          PropertiesAttribute(Map("foo" -> "bar")),
          DetailsAttribute(Map("whatsNew" -> VString("v1.0: First release"))),
          TitleAttribute("Workflow with metadata"),
          TypesAttribute(Vector("calculator")),
          SummaryAttribute("A workflow that defines some metadata")
      )
  }

  it should "recognize workflow parameter metadata" in {
    val path = pathFromBasename("compiler", "wf_param_meta.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case _                  => throw new Exception("primaryCallable is not a workflow")
    }
    val input_cvars: Vector[Parameter] = workflow.inputs.map {
      case (c: Parameter, _) => c
      case _                 => throw new Exception("Invalid workflow input ${other}")
    }
    input_cvars.sortWith(_.name < _.name) shouldBe Vector(
        Parameter(
            WdlDxName.fromSourceName("x"),
            TInt,
            Some(VInt(3)),
            Vector(
                LabelAttribute("Left-hand side"),
                DefaultAttribute(VInt(3))
            )
        ),
        Parameter(
            WdlDxName.fromSourceName("y"),
            TInt,
            Some(VInt(5)),
            Vector(
                LabelAttribute("Right-hand side"),
                DefaultAttribute(VInt(5))
            )
        )
    )
  }

  it should "handle adjunct files in workflows and tasks" in {
    val path = pathFromBasename("compiler", "wf_readme.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]

    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val workflow = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case _                  => throw new Exception("primaryCallable is not a workflow")
    }
    workflow.attributes match {
      case Vector(DescriptionAttribute(desc)) =>
        desc shouldBe "This is the readme for the wf_linear workflow."
      case _ =>
        throw new Exception("Expected one workflow description attribute")
    }

    val addApp = getApplicationByName("add", bundle)
    addApp.attributes.size shouldBe 2
    addApp.attributes.foreach {
      case DescriptionAttribute(text) =>
        text shouldBe "This is the readme for the wf_linear add task."
      case DeveloperNotesAttribute(text) =>
        text shouldBe "Developer notes defined in WDL"
      case other => throw new Exception(s"Invalid  for add task ${other}")
    }

    val mulApp = getApplicationByName("mul", bundle)
    mulApp.attributes match {
      case Vector(DescriptionAttribute(text)) =>
        text shouldBe "Description defined in WDL"
      case other =>
        throw new Exception(s"expected one description attribute, not ${other}")
    }

    val incApp = getApplicationByName("inc", bundle)
    incApp.attributes.size shouldBe 0
  }

  it should "work correctly with pairs as call inputs in a scatter" taggedAs EdgeTest in {
    val path = pathFromBasename("subworkflows", basename = "scatter_subworkflow_with_optional.wdl")
    val cFlagsNotQuiet = cFlags.filter(_ != "-quiet")
    val args = path.toString :: cFlagsNotQuiet
    //          "--verbose" ::
    //          :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "work correctly with pairs in a simple scatter" taggedAs EdgeTest in {
    val path = pathFromBasename("frag_runner", basename = "scatter_with_eval.wdl")
    val cFlagsNotQuiet = cFlags.filter(_ != "-quiet")
    val args = path.toString :: cFlagsNotQuiet
    //          :: "--verbose"
    //          :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "work correctly with nested pairs in a simple scatter" taggedAs EdgeTest in {
    val path = pathFromBasename("bugs", basename = "apps-370.wdl")
    val cFlagsNotQuiet = cFlags.filter(_ != "-quiet")
    val args = path.toString :: cFlagsNotQuiet
    //          :: "--verbose"
    //          :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  it should "work correctly with a complex scatter" taggedAs EdgeTest in {
    val path = pathFromBasename("bugs", basename = "apps-378.wdl")
    val cFlagsNotQuiet = cFlags.filter(_ != "-quiet")
    val args = path.toString :: cFlagsNotQuiet
    //          :: "--verbose"
    //          :: "--verboseKey" :: "GenerateIR"
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
  }

  // Check parameter_meta pattern: ["array"]
  it should "recognize pattern in parameters_meta via WDL" in {
    val path = pathFromBasename("compiler", "pattern_params.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepTask = getApplicationByName("pattern_params_cgrep", bundle)
    cgrepTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "in_file" -> Vector(
            HelpAttribute("The input file to be searched"),
            PatternsArray(Vector("*.txt", "*.tsv")),
            GroupAttribute("Common"),
            LabelAttribute("Input file")
        ),
        "pattern" -> Vector(
            HelpAttribute("The pattern to use to search in_file"),
            GroupAttribute("Common"),
            LabelAttribute("Search pattern")
        ),
        "out_file" -> Vector(
            PatternsArray(
                Vector("*.txt", "*.tsv")
            ),
            GroupAttribute("Common"),
            LabelAttribute("Output file")
        )
    )
  }

  // Check parameter_meta pattern: {"object"}
  it should "recognize pattern object in parameters_meta via WDL" in {
    val path = pathFromBasename("compiler", "pattern_obj_params.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepTask = getApplicationByName("pattern_params_obj_cgrep", bundle)
    cgrepTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "in_file" -> Vector(
            HelpAttribute("The input file to be searched"),
            PatternsAttribute(
                PatternsObject(
                    Vector("*.txt", "*.tsv"),
                    Some("file"),
                    Vector("foo", "bar")
                )
            ),
            GroupAttribute("Common"),
            LabelAttribute("Input file")
        ),
        "pattern" -> Vector(
            HelpAttribute("The pattern to use to search in_file"),
            GroupAttribute("Common"),
            LabelAttribute("Search pattern")
        ),
        "out_file" -> Vector(
            PatternsArray(Vector("*.txt", "*.tsv")),
            GroupAttribute("Common"),
            LabelAttribute("Output file")
        )
    )
  }

  it should "recognize help, group, and label in parameters_meta via WDL" in {
    val path = pathFromBasename("compiler", "help_input_params.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }

    val cgrepTask = getApplicationByName("help_input_params_cgrep", bundle)
    cgrepTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "in_file" -> Vector(
            HelpAttribute("The input file to be searched"),
            GroupAttribute("Common"),
            LabelAttribute("Input file")
        ),
        "pattern" -> Vector(
            DescriptionAttribute("The pattern to use to search in_file"),
            GroupAttribute("Common"),
            LabelAttribute("Search pattern")
        ),
        "s" -> Vector(HelpAttribute("This is help for s"))
    )

    val diffTask = getApplicationByName("help_input_params_diff", bundle)
    diffTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "a" -> Vector(
            HelpAttribute("lefthand file"),
            GroupAttribute("Files"),
            LabelAttribute("File A")
        ),
        "b" -> Vector(
            HelpAttribute("righthand file"),
            GroupAttribute("Files"),
            LabelAttribute("File B")
        )
    )
  }

  // `paramter: "stream"` should not be converted to an attribute - it is only accessed at runtime
  it should "ignore stream attribute" in {
    val path = pathFromBasename("compiler", "streaming_files.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val cgrepTask = getApplicationByName("cgrep", bundle)
    cgrepTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "in_file" -> Vector()
    )
    val diffTask = getApplicationByName("diff", bundle)
    diffTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "a" -> Vector(),
        "b" -> Vector()
    )
  }

  // `parameter: {stream: true}` should not be converted to an attribute - it is only accessed at runtime
  it should "ignore the streaming object annotation" in {
    val path = pathFromBasename("compiler", "streaming_files_obj.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val cgrepTask = getApplicationByName("cgrep", bundle)
    cgrepTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "in_file" -> Vector()
    )
    val diffTask = getApplicationByName("diff", bundle)
    diffTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "a" -> Vector(),
        "b" -> Vector()
    )
  }

  it should "ignore the streaming annotation for wdl draft2" in {
    val path = pathFromBasename("draft2", "streaming.wdl")
    val args = path.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    val diffTask = getApplicationByName("diff", bundle)
    diffTask.inputs.map(param => param.name -> param.attributes).iterator sameElements Vector(
        "a" -> Vector(),
        "b" -> Vector()
    )
  }

  it should "recognize default and per-workflow attributes" in {
    val path = pathFromBasename("nested", "four_levels.wdl")
    val extraPath = pathFromBasename("nested/extras", "four_levels_extras1.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("unexpected")
    }
    // there are two scatters - one should have its chunkSize set
    // by the global default and the other should have it set by
    // a per-scatter override
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case _                  => throw new Exception("missing primary workflow")
    }
    val scatter = wf.stages(1)
    val outerScatterApplet = bundle.allCallables.get(scatter.calleeName) match {
      case Some(applet: Application) => applet
      case _                         => throw new Exception(s"missing applet for stage ${scatter}")
    }
    val innerName = outerScatterApplet.kind match {
      case frag: ExecutableKindWfFragment =>
        frag.scatterChunkSize shouldBe Some(10)
        frag.call match {
          case Some(innerName) => innerName
          case _               => throw new Exception(s"wrong calls ${frag.call}")
        }
      case _ => throw new Exception(s"wrong applet kind ${outerScatterApplet.kind}")
    }
    val innerScatterApplet = bundle.allCallables.get(innerName) match {
      case Some(applet: Application) => applet
      case _                         => throw new Exception(s"missing applet for stage ${innerName}")
    }
    innerScatterApplet.kind match {
      case frag: ExecutableKindWfFragment =>
        frag.scatterChunkSize shouldBe Some(3)
      case _ => throw new Exception(s"wrong applet kind ${innerScatterApplet.kind}")
    }
  }

  it should "recognize per-workflow attributes with global default" in {
    val path = pathFromBasename("nested", "four_levels.wdl")
    val extraPath = pathFromBasename("nested/extras", "four_levels_extras2.json")
    val args = path.toString :: "--extras" :: extraPath.toString :: cFlags
    val retval = Main.compile(args.toVector)
    retval shouldBe a[SuccessfulCompileIR]
    val bundle = retval match {
      case SuccessfulCompileIR(ir) => ir
      case _                       => throw new Exception("wrong termination type")
    }
    // there are two scatters - one should have its chunkSize set
    // by the global default and the other should have it set by
    // a per-scatter override
    val stages = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf.stages
      case _                  => throw new Exception("missing primary workflow")
    }
    val scatter = stages(1)
    val outerScatterApplet = bundle.allCallables.get(scatter.calleeName) match {
      case Some(applet: Application) => applet
      case _                         => throw new Exception(s"missing applet for stage ${scatter}")
    }
    val innerName = outerScatterApplet.kind match {
      case frag: ExecutableKindWfFragment =>
        frag.scatterChunkSize shouldBe Some(Constants.JobPerScatterDefault)
        frag.call match {
          case Some(innerName) => innerName
          case _               => throw new Exception(s"wrong calls ${frag.call}")
        }
      case _ => throw new Exception(s"wrong applet kind ${outerScatterApplet.kind}")
    }
    val innerScatterApplet = bundle.allCallables.get(innerName) match {
      case Some(applet: Application) => applet
      case _                         => throw new Exception(s"missing applet for stage ${innerName}")
    }
    innerScatterApplet.kind match {
      case frag: ExecutableKindWfFragment =>
        frag.scatterChunkSize shouldBe Some(3)
      case _ => throw new Exception(s"wrong applet kind ${innerScatterApplet.kind}")
    }
  }

  it should "translate a simple CWL tool" in {
    val path = pathFromBasename("cwl", "cat.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) =>
        bundle.allCallables.size shouldBe 1
        val applet = bundle.primaryCallable match {
          case Some(applet: Application) => applet
          case other =>
            throw new AssertionError(s"expected primaryCallable to be an applet, not ${other}")
        }
        applet.name shouldBe "cat"
        applet.kind shouldBe ExecutableKindApplet
        applet.attributes shouldBe Vector(DescriptionAttribute("Write a file to stdout using cat"))
        applet.container shouldBe NoImage
        applet.inputs shouldBe Vector(Parameter(CwlDxName.fromSourceName("file"), TFile),
                                      TargetParam)
        applet.outputs shouldBe Vector(Parameter(CwlDxName.fromSourceName("contents"), TFile))
      case other =>
        throw new AssertionError(s"expected SuccessfulCompileIR, not ${other}")
    }
  }

  it should "translate a simple CWL tool with JS expressions" in {
    val path = pathFromBasename("cwl", "params.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a simple CWL tool with Any type input" in {
    val path = pathFromBasename("cwl", "params2.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a CWL tool with inline record type" in {
    val path = pathFromBasename("cwl", "record-in-format.cwl.json")
    val args = path.toString :: cFlags
    val result = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        throw new AssertionError(s"expected SuccessfulCompileIR, not ${other}")
    }
    val tool = result.primaryCallable match {
      case Some(tool: Application) => tool
      case other =>
        throw new Exception(s"expected compilation result to be a tool, not ${other}")
    }
    tool.inputs
      .collectFirst {
        case inp if inp.name.decoded == "record_input" => inp
      }
      .getOrElse(throw new Exception("missing field record_input"))
  }

  it should "translate a workflow with scatter inside conditional" in {
    val path = pathFromBasename("bugs", "scatter_inside_if.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a simple CWL workflow" in {
    val path = pathFromBasename("cwl", "count-lines1-wf.cwl.json")
    val args = path.toString :: cFlags
    val result = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        throw new Exception(s"expected success not ${other}")
    }
    val wf = result.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other =>
        throw new Exception(s"expected Workflow not ${other}")
    }
    wf.inputs shouldBe Vector(
        Parameter(CwlDxName.fromSourceName("file1"), TFile) -> StageInputWorkflowLink(
            Parameter(CwlDxName.fromSourceName("file1"), TFile)
        )
    )
    wf.outputs shouldBe Vector(
        Parameter(CwlDxName.fromSourceName("count_output"), TInt) -> StageInputStageLink(
            DxWorkflowStage("stage-1"),
            Parameter(CwlDxName.fromSourceName("output"), TInt)
        )
    )
    wf.stages.size shouldBe 2
    wf.stages(0).dxStage.id shouldBe "stage-0"
    wf.stages(0).calleeName shouldBe "wc-tool"
    wf.stages(0).inputs shouldBe Vector(
        StageInputWorkflowLink(Parameter(CwlDxName.fromSourceName("file1"), TFile)),
        StageInputStatic(VString("step1"))
    )
    wf.stages(0).outputs shouldBe Vector(
        Parameter(CwlDxName.fromSourceName("output"), TFile)
    )
    wf.stages(1).dxStage.id shouldBe "stage-1"
    wf.stages(1).calleeName shouldBe "parseInt-tool"
    wf.stages(1).inputs shouldBe Vector(
        StageInputStageLink(DxWorkflowStage("stage-0"),
                            Parameter(CwlDxName.fromSourceName("output"), TFile)),
        StageInputStatic(VString("step2"))
    )
    wf.stages(1).outputs shouldBe Vector(
        Parameter(CwlDxName.fromSourceName("output"), TInt)
    )
  }

  it should "translate a packed CWL workflow" in {
    val path = pathFromBasename("cwl", "any-type-compat.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        throw new Exception(s"expected success not ${other}")
    }
  }

  it should "translate a CWL workflow with nested workflow" in {
    val path = pathFromBasename("cwl", "count-lines18-wf.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other =>
        throw new Exception(s"expected success not ${other}")
    }
  }

  it should "translate a workflow with struct field as call argument" in {
    val path = pathFromBasename("struct", "struct_field.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a workflow with references between output variables" in {
    val path = pathFromBasename("draft2", "output_references.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a CWL workflow with JavaScript expressions" in {
    val path = pathFromBasename("cwl", "timelimit2-wf.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a CWL workflow with multiple scatter sources" in {
    val path = pathFromBasename("cwl", "scatter-valuefrom-wf2.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a CWL workflow with auto-generated embedded process IDs" in {
    val path = pathFromBasename("cwl", "cond-wf-003_nojs.cwl.json")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => ()
      case other                       => throw new Exception(s"expected success not ${other}")
    }
  }

  it should "translate a CWL workflow with workflow input link of different name" in {
    val path = pathFromBasename("cwl", "cond-wf-011_nojs.cwl.json")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other                       => throw new Exception(s"expected success not ${other}")
    }
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other              => throw new Exception(s"expected workflow not ${other}")
    }
    val stageApplet = bundle.allCallables(wf.stages.head.calleeName) match {
      case applet: Application => applet
      case other               => throw new Exception(s"expected applet not ${other}")
    }
    val stageParams = stageApplet.inputs.map(param => param.name -> param).toMap
    stageParams.keySet
      .map(_.decodedIdentifier) shouldBe Set("in1", "in2", "in3", "test")
    val targetApplet = stageApplet.kind match {
      case ExecutableKindWfFragment(Some(call), _, _, _) => bundle.allCallables(call)
      case other =>
        throw new Exception(s"expetected fragment not ${other}")
    }
    val targetParams = targetApplet.inputVars.map(param => param.name -> param).toMap
    targetParams.keySet
      .map(_.decodedIdentifier) shouldBe Set("in1", "in2", "in3", "target")
  }

  it should "translate a CWL workflow with stage input with linkMerge" in {
    val path = pathFromBasename("cwl", "count-lines19-wf.cwl.json")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other                       => throw new Exception(s"expected success not ${other}")
    }
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other              => throw new Exception(s"expected workflow not ${other}")
    }
    val stageApplet = bundle.allCallables(wf.stages.head.calleeName) match {
      case applet: Application => applet
      case other               => throw new Exception(s"expected applet not ${other}")
    }
    stageApplet.kind should matchPattern {
      case ExecutableKindWfFragment(Some("wc3-tool"), Vector(0), _, None) =>
    }
  }

  it should "translate a WDL workflow with dx runtime attributes" in {
    val path = pathFromBasename("bugs", "dx_runtime_keys.wdl")
    val args = path.toString :: cFlags
    Main.compile(args.toVector) shouldBe a[SuccessfulCompileIR]
  }

  it should "translate a CWL workflow with a multi-type default value" in {
    val path = pathFromBasename("cwl", "io-union-input-default-wf.cwl.json")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other                       => throw new Exception(s"expected success not ${other}")
    }
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other              => throw new Exception(s"expected workflow not ${other}")
    }
    wf.inputs.size shouldBe 1
    wf.inputs.head._1.defaultValue shouldBe Some(VString("the default value"))
  }

  it should "translate a CWL workflow with a step input source and default value" in {
    val path = pathFromBasename("cwl", "dynresreq-workflow-stepdefault.cwl.json")
    val args = path.toString :: cFlags
    val bundle = Main.compile(args.toVector) match {
      case SuccessfulCompileIR(bundle) => bundle
      case other                       => throw new Exception(s"expected success not ${other}")
    }
    val wf = bundle.primaryCallable match {
      case Some(wf: Workflow) => wf
      case other              => throw new Exception(s"expected workflow not ${other}")
    }
    wf.stages.size shouldBe 2
    val stage0Applet = bundle.allCallables(wf.stages(0).calleeName) match {
      case applet: Application => applet
      case other               => throw new Exception(s"expected applet not ${other}")
    }
    stage0Applet.kind should matchPattern {
      case ExecutableKindWfFragment(_, _, _, _) =>
    }
  }
}
