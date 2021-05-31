package dx.compiler

import dx.core.ir.{Parameter, Value}
import dx.translator.ParameterAttributes

package object CompilerUtils {

  // Return a parameter's dependencies on non-local files
  def fileDependenciesFromParam(param: Parameter): Set[String] = {
    val default = param.defaultValue match {
      case Some(Value.VFile(str)) => Vector(str)
      case _                      => Vector.empty
    }
    val attrs = param.attributes
      .flatMap(attr =>
        attr match {
          case ParameterAttributes.DefaultAttribute(Value.VFile(str)) => Vector(str)
          case ParameterAttributes.ChoicesAttribute(choices) =>
            choices
              .collect {
                case ParameterAttributes.SimpleChoice(Value.VFile(str)) => str
                case ParameterAttributes.FileChoice(str, _)             => str
              }
          case ParameterAttributes.SuggestionsAttribute(suggestions) =>
            suggestions
              .collect {
                case ParameterAttributes.SimpleSuggestion(Value.VFile(str)) => str
                case ParameterAttributes.FileSuggestion(Some(str), _, _, _) => str
              }
          case _ => Vector.empty
        }
      )
      .toVector

    (default ++ attrs).filter(s => s.contains("://")).toSet
  }
}
