import "http://localhost:8000/private" as subworkflow

workflow public_http_import {

  Array[Int] ts = range(3)
  call subworkflow.subwf as subwfT { input: is = ts }

  Array[Int] fs = subwfT.js
  call subworkflow.subwf as subwfF { input: is = fs }

  output {
    Array[Int] initial = ts
    Array[Int] intermediate = subwfT.js
    Array[Int] result = subwfF.js
  }
}
