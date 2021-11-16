import "subworkflow_wt.wdl" as subworkflow

workflow aliased_subworkflows {

  Array[Int] ts = range(3)
  call subworkflow.subwf_wt as subwfT { input: is = ts }

  Array[Int] fs = subwfT.js
  call subworkflow.subwf_wt as subwfF { input: is = fs }

  output {
    Array[Int] initial = ts
    Array[Int] intermediate = subwfT.js
    Array[Int] result = subwfF.js
  }
}
