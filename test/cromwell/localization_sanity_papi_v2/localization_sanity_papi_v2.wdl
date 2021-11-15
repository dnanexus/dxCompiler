version 1.0

task make_files {
  command <<<
    names=(a b c)
    mkdir -p "${names[@]}"
    for name in "${names[@]}"; do
      touch "${name}/dummy.txt" # the first file is not bulk-transferred via `gsutil cp -I ...` which is what this test is about.
      touch "${name}/${name}.txt"
    done
  >>>
  output {
    # Intentionally not globbed as the current implementaton of globbing would defeat what this test
    # is trying to assert.
    Array[File] files = ["a/dummy.txt", "a/a.txt", "b/dummy.txt", "b/b.txt", "c/dummy.txt", "c/c.txt"]
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

task sanity_check {
  input {
    Array[File] files
  }
  command <<<
    names=(a b c)
    for name in "${names[@]}"; do
      file="${name}.txt"
      echo "file $file: $(find . -name $file | wc -l)"
    done
  >>>
  output {
    Array[String] lines = read_lines(stdout())
  }
  runtime {
    docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
  }
}

workflow localization_sanity {
  call make_files
  call sanity_check { input: files = make_files.files }
}
