package dx.core

import java.nio.file.{Files, Path, Paths}

import dx.util.{Logger, TraceLevel}

import scala.reflect.ClassTag

object CliUtils {
  case class OptionParseException(message: String) extends Exception(message)

  // simple command line parsing - only handles three types of options:
  // * flag
  // * one or more values, specified once
  // * single value, specified any number of times
  sealed trait Opt
  case class Flag(enabled: Boolean = true) extends Opt
  case class SingleValueOption[T](value: T) extends Opt
  case class ListOption[T](value: Vector[T]) extends Opt

  case class Options(options: Map[String, Opt]) {
    def update(name: String, value: Opt): Options = {
      Options(options + (name -> value))
    }

    def contains(name: String): Boolean = options.contains(name)

    def get(name: String): Option[Opt] = options.get(name)

    def getFlag(name: String, default: Boolean = false): Boolean = {
      options.get(name) match {
        case Some(Flag(b)) => b
        case None          => default
        case _ =>
          throw OptionParseException(s"Option ${name} is not a flag")
      }
    }

    def getValue[T: ClassTag](name: String): Option[T] = {
      options.get(name).map {
        case SingleValueOption(value: T) => value
        case other =>
          throw OptionParseException(s"Unexpected value ${other} to option ${name}")
      }
    }

    def getRequiredValue[T: ClassTag](name: String): T = {
      getValue[T](name).getOrElse(
          throw OptionParseException(s"Missing required option ${name}")
      )
    }

    def getValueOrElse[T: ClassTag](name: String, default: T): T = {
      getValue[T](name).getOrElse(default)
    }

    // TODO: is there a way to do this without resorting to asInstanceOf?
    def getList[T: ClassTag](name: String): Vector[T] = {
      options.get(name) match {
        case Some(list: ListOption[_]) => list.value.asInstanceOf[Vector[T]]
        case None                      => Vector.empty[T]
        case other =>
          throw OptionParseException(s"Unexpected value ${other} to option ${name}")
      }
    }
  }

  object Options {
    lazy val empty: Options = Options(Map.empty[String, Opt])
  }

  /**
    * Specification for a command line option to parse.
    */
  trait OptionSpec {

    /**
      * Whether the option may be specified multiple times
      */
    val multiple: Boolean = false

    /**
      * Alias for this option
      */
    val alias: Option[String] = None

    /**
      * Parse the option and return the real name (e.g. if this
      * option is aliased) and value.
      * @param name the user-specified option name
      * @param values the command-line values
      * @return
      */
    def parse(name: String, values: Vector[String], curValue: Option[Opt]): (String, Opt) = {
      if (curValue.isDefined && !multiple) {
        throw OptionParseException(s"Option ${name} appears multiple times")
      }
      (alias.getOrElse(name), parseValues(name, values, curValue))
    }

    def parseValues(name: String, values: Vector[String], curValue: Option[Opt]): Opt
  }

  /**
    * Specification for a flag (boolean) option
    */
  case class FlagOptionSpec(override val alias: Option[String] = None) extends OptionSpec {
    def parseValues(name: String, values: Vector[String], curValue: Option[Opt]): Flag = {
      if (curValue.nonEmpty) {
        throw OptionParseException(s"Option ${name} specified multiple times")
      }
      if (values.nonEmpty) {
        throw OptionParseException(s"Flag option ${name} has unexpected values ${values}")
      }
      Flag()
    }
  }

  object FlagOptionSpec {
    val default: FlagOptionSpec = FlagOptionSpec()
  }

  abstract class SingleValueOptionSpec[T](override val alias: Option[String] = None,
                                          choices: Vector[T] = Vector.empty,
                                          override val multiple: Boolean = false)
      extends OptionSpec {
    def parseValues(name: String, values: Vector[String], curValue: Option[Opt]): Opt = {
      if (values.size != 1) {
        throw OptionParseException(
            s"Expected option ${name} to have 1 value, found ${values}"
        )
      }
      if (curValue.nonEmpty && !multiple) {
        throw OptionParseException(s"Option ${name} specified multiple times")
      }
      val value = parseValue(values.head)
      if (choices.nonEmpty && !choices.contains(value)) {
        throw OptionParseException(
            s"Unexpected value ${value} to option ${name}"
        )
      }
      (curValue, multiple) match {
        case (Some(ListOption(values: Vector[T])), true) =>
          ListOption[T](values :+ value)
        case (None, true) =>
          ListOption[T](Vector(value))
        case (None, false) =>
          SingleValueOption[T](value)
        case _ =>
          throw OptionParseException(s"Unexpected value ${value} to option ${name}")
      }
    }

    def parseValue(value: String): T
  }

  /**
    * Specification for a String option that takes a single value
    * @param choices (optional) set of allowed values
    */
  case class StringOptionSpec(override val alias: Option[String] = None,
                              choices: Vector[String] = Vector.empty,
                              override val multiple: Boolean = false)
      extends SingleValueOptionSpec[String](alias, choices, multiple) {
    def parseValue(value: String): String = value
  }

  object StringOptionSpec {
    lazy val one: StringOptionSpec = StringOptionSpec()
    lazy val list: StringOptionSpec = StringOptionSpec(multiple = true)
  }

  case class IntOptionSpec(override val alias: Option[String] = None,
                           choices: Vector[Int] = Vector.empty,
                           override val multiple: Boolean = false)
      extends SingleValueOptionSpec[Int](alias, choices, multiple) {
    def parseValue(value: String): Int = value.toInt
  }

  object IntOptionSpec {
    lazy val one: IntOptionSpec = IntOptionSpec()
    lazy val list: IntOptionSpec = IntOptionSpec(multiple = true)
  }

  case class PathOptionSpec(mustExist: Boolean = false,
                            override val alias: Option[String] = None,
                            override val multiple: Boolean = false)
      extends SingleValueOptionSpec[Path](alias, multiple = multiple) {
    def parseValue(value: String): Path = {
      val path = Paths.get(value)
      if (Files.exists(path)) {
        path.toAbsolutePath
      } else if (mustExist) {
        throw OptionParseException(s"Path ${path} does not exist")
      } else {
        path
      }
    }
  }

  object PathOptionSpec {
    lazy val default: PathOptionSpec = PathOptionSpec()
    lazy val mustExist: PathOptionSpec = PathOptionSpec(mustExist = true)
    lazy val listMustExist: PathOptionSpec = PathOptionSpec(mustExist = true, multiple = true)
  }

  type InternalOptions = Map[String, OptionSpec]

  val QuietOpt = "quiet"
  val VerboseOpt = "verbose"
  val VerboseKeyOpt = "verboseKey"
  val TraceLevelOpt = "traceLevel"
  val LogFileOpt = "logFile"

  private val SimpleOptions: InternalOptions = Map(
      QuietOpt -> FlagOptionSpec.default,
      VerboseOpt -> FlagOptionSpec.default,
      VerboseKeyOpt -> StringOptionSpec.list,
      TraceLevelOpt -> IntOptionSpec.one.copy(choices = Vector(0, 1, 2)),
      LogFileOpt -> PathOptionSpec.default
  )

  // Split arguments into sub-lists, one per each option.
  // For example:
  //    --sort relaxed --reorg --compile-mode IR
  // =>
  //    [[--sort, relaxed], [--reorg], [--compile-mode, IR]]
  //
  def parseCommandLine(args: Vector[String],
                       specs: InternalOptions,
                       deprecated: Set[String] = Set.empty): Options = {
    def isOption(word: String): Boolean = word.startsWith("-")

    val OptionRegexp = "^-+(.+)$".r

    def normOption(word: String): String = {
      // normalize a keyword, remove leading dashes
      // letters to lowercase.
      //
      // "--Archive" -> "archive"
      // "--archive-only -> "archiveonly"
      word match {
        case OptionRegexp(rest) => rest
        case _ =>
          throw OptionParseException(s"${word} is not an option")
      }
    }

    val allSpecs = SimpleOptions ++ specs

    def createOpt(args: Vector[String], curOpts: Options): Options = {
      val name = args.headOption.getOrElse(throw new Exception("Option with no name"))
      if (deprecated.contains(name)) {
        Logger.warning(s"Option ${name} is deprecated and will be ignored")
        return curOpts
      }
      val values = args.tail
      val curValue = curOpts.get(name)
      val spec = allSpecs.getOrElse(
          name,
          throw new Exception(s"Unexpected option ${name}")
      )
      val (realName, opt) = spec.parse(name, values, curValue)
      curOpts.update(realName, opt)
    }

    val (lastOpt, optMap) = args.foldLeft(Vector.empty[String], Options.empty) {
      case ((curOpt, allOpts), word) if isOption(word) && curOpt.nonEmpty =>
        val newOpts = createOpt(curOpt, allOpts)
        (Vector(normOption(word)), newOpts)
      case ((_, allOpts), word) if isOption(word) =>
        (Vector(normOption(word)), allOpts)
      case ((curOpt, allOpts), word) if curOpt.nonEmpty =>
        (curOpt :+ word, allOpts)
      case _ =>
        throw new Exception("Option must precede value")
    }

    if (lastOpt.nonEmpty) {
      createOpt(lastOpt, optMap)
    } else {
      optMap
    }
  }

  // logging

  def initLogger(options: Options): Logger = {
    val verboseKeys: Set[String] = options.getList[String](VerboseKeyOpt).toSet
    val traceLevel = options
      .getValue[Int](TraceLevelOpt)
      .getOrElse(
          if (options.getFlag(VerboseOpt)) {
            TraceLevel.Verbose
          } else {
            TraceLevel.None
          }
      )
    val logFile = options.getValue[Path]("logFile")
    val logger = Logger(quiet = options.getFlag(QuietOpt),
                        traceLevel = traceLevel,
                        keywords = verboseKeys,
                        logFile = logFile)
    Logger.set(logger)
    logger
  }

  // command exit states

  def exceptionToString(e: Throwable): String = {
    val sw = new java.io.StringWriter
    e.printStackTrace(new java.io.PrintWriter(sw))
    sw.toString
  }

  private def failureMessage(message: String, exception: Option[Throwable]): String = {
    val exStr = exception.map(exceptionToString)
    (message, exStr) match {
      case (s, Some(e)) if s.nonEmpty =>
        s"""${message}
           |${e}""".stripMargin
      case ("", Some(e)) => e
      case (s, None)     => s
    }
  }

  trait Termination {
    def message: String
  }

  trait SuccessfulTermination extends Termination
  case class Success(override val message: String = "") extends SuccessfulTermination

  trait UnsuccessfulTermination extends Termination {
    val message: String
    val exception: Option[Throwable] = None
    override def toString: String = failureMessage(message, exception)
  }
  case class Failure(override val message: String = "",
                     override val exception: Option[Throwable] = None)
      extends UnsuccessfulTermination
  case class BadUsageTermination(override val message: String = "",
                                 override val exception: Option[Throwable] = None)
      extends UnsuccessfulTermination

  def terminate(termination: Termination, usageMessage: String): Unit = {
    val rc = termination match {
      case _: SuccessfulTermination =>
        if (termination.message.nonEmpty) {
          println(termination.message)
        }
        0
      case BadUsageTermination("", e) =>
        Logger.error(failureMessage(usageMessage, e))
        1
      case term: UnsuccessfulTermination =>
        val msg = failureMessage(term.message, term.exception)
        if (msg.nonEmpty) {
          Logger.error(msg)
        }
        1
    }
    System.exit(rc)
  }
}
