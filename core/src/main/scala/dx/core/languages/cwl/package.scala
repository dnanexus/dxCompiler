package dx.core.languages

import dx.core.Constants
import dx.core.ir.{DxName, Parameter, SimpleDxName}
import dx.core.ir.Type.{TOptional, TString}

package object cwl {
  val Target: DxName = SimpleDxName.fromSourceName("target", Some(Constants.ComplexValueKey))
  val TargetParam: Parameter = Parameter(Target, TOptional(TString))
}
