package dx.core.languages

import dx.core.Constants
import dx.core.ir.{DxName, Parameter}
import dx.core.ir.Type.{TOptional, TString}

package object cwl {
  val Target: DxName = CwlDxName.fromSourceName("target", None, Some(Constants.ComplexValueKey))
  val TargetParam: Parameter = Parameter(Target, TOptional(TString))
}
