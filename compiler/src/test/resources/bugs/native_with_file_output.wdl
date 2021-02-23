version 1.0

import "subwf.wdl" as sub

workflow native_with_file_output {
  input {
    Array[String]+ f_urls
    String target_s3
    File config_file
    String up_dir
  }

  call aws_s3_to_platform_files {
    input: f_urls = f_urls, target_s3 = target_s3, config_file = config_file, up_dir = up_dir
  }

  call sub.foo {
    input: inp = aws_s3_to_platform_files.transferred_files
  }

  output {
    Array[File]+ transferred_files = foo.out
  }
}

task aws_s3_to_platform_files {
  input {
    Array[String]+ f_urls
    String target_s3
    File config_file
    String up_dir
    Int worker_max = 1
    String bandwidth = "alot"
    String? additional_upload_param
    Boolean upload_direct_to_proj = true
  }

  command <<< >>>

  output {
    File upload_report = "placeholder.txt"
    Array[File]+ transferred_files = ["placeholder.txt"]
  }

  meta {
    type: "native"
    id: "app-F0px8Y005x8p3Z3B8pzx3VfZ"
  }
}
