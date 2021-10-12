 version 1.1
 workflow test_filter {
  input { 
  Array[String]+ sample_names = ["A", "B", "C"] 
  Array[Int]+ numbers = [1, 2, 3] 
  Array[String] samples_to_exclude } 
  scatter (sample_int in zip(sample_names, numbers)) { 
      scatter (exclude in samples_to_exclude) {
        Boolean? drop = if sample_int.left == exclude then true else None
      }
      Pair[String, Int]? kept_sample = if length(select_all(drop)) == 0 then sample_int else None
  } 
  Array[Pair[String, Int]] kept_samples = select_all(kept_sample) Pair[Array[String], Array[Int]] samples_ints = unzip(kept_samples) 
  output { 
  Array[String] samples_passing = samples_ints.left
 Array[Int] ints_passing = samples_ints.right 
 }
 } 