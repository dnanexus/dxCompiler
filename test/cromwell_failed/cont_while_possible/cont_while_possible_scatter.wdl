task odd_numbers_only {
Int i
  command {
	  if [ $((${i}%2)) -eq 0 ]; then
	    exit 1;
	  fi
	}
	runtime {
		docker: "dx://file-G66qpGj0yzZq02K9313pJg5G"
	}
}

workflow cont_while_possible_scatter {
  Array[Int] myRange = range(4)

	scatter(i in myRange){
	  call odd_numbers_only {input: i = i}
	}
}
