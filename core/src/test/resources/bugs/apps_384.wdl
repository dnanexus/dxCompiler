version 1.0

workflow apps_384 {
  input {
    String dataset = "dataset"
    Array[Pair[String, Pair[String, String]]] manifest = [("str1", ("str2", "str3"))]
    Pair[String, Pair[String, String]] row = manifest[0]
  }

  String RG_LB = "${dataset}_${dataset}_${row.left}"

  call cat {
    input:
      rg = "'@RG\\tID:" + RG_LB + "\\tPL:Illumina\\tPU:" + RG_LB + "\\tLB:" + RG_LB + "\\tSM:" + row.left + "'"
  }

  output {
    String dataset_out = dataset
    String row1 = manifest[0].left
    String row2 = manifest[0].right.left
    String row3 = manifest[0].right.right
    String rg = cat.rg_out
  }
}

task cat {
  input {
    String rg
  }

  command {}

  output {
    String rg_out = rg
  }
}
