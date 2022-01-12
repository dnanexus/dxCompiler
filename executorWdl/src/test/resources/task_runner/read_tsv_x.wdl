# Generate sets of intervals for scatter-gathering over chromosomes
task read_tsv_x {
    # Use python to create a string parsed into a wdl Array[Array[String]]
    command<<<
    python <<CODE
    tsv_list = []
    ll = [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    for l in ll:
      tsv_list.append('\t'.join(l))
    tsv_string = '\n'.join(tsv_list)
    print(tsv_string)
    CODE
    >>>

    output {
      Array[Array[String]] result = read_tsv(stdout())
    }
}
