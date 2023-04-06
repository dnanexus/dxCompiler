# A file with one task should compile to one runnable
# applet.

task spaces_in_file_paths {
    command {
        echo "A" > 'file A.txt'
        echo "B" > 'file B.txt'
        echo "C" > 'file C.txt'
    }
    output {
        Array[File] texts = glob("*.txt")
    }
}
