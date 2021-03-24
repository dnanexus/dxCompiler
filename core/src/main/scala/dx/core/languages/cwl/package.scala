package dx.core.languages

import dx.core.ir.Parameter
import dx.core.ir.Type.{TOptional, TString}

package object cwl {
  val Target = "target___"
  val TargetParam: Parameter = Parameter(Target, TOptional(TString))
}
