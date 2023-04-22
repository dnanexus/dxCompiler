import "../imports/subworkflow_ii.wdl" as subworkflow_ii

workflow subworkflows_in_ifs {

  if (true) {
    Array[Int] ts = range(3)
    call subworkflow_ii.subwf_ii as subwfTrue { input: is = ts }
  }
  if (false) {
    Array[Int] fs = range(3)
    call subworkflow_ii.subwf_ii as subwfFalse { input: is = fs }
  }

  output {
    Array[Int]? tjs = subwfTrue.js
    Array[Int]? fjs = subwfFalse.js
  }
}
