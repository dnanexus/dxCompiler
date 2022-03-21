version 1.1

task reuse_multiply {
	input {
		Int first
		Int second
	}

	command <<< >>>

	runtime {
		container: "ubuntu:20.04"
	}

	Int product = first * second

	output {
		Int out = product
	}
}

task reuse_print {
	input {
		String str
	}

	command <<<
	set -euxo pipefail
	echo ~{str} > reuse_print.txt
	>>>

	runtime {
		container: "ubuntu:20.04"
	}

	output {
		File out = "reuse_print.txt"
	}
}

workflow reuse {

	# Frag applet via sub-expression
	# declaration of inputs is different from V1. Therefore this frag should be rebuilt
	Int f = 1 + 10000
	Int s = 4 + 5
	call reuse_multiply as sub_expr {
		input: first = f, second = s
	}

	# Frag applet via conditional, single call
	Boolean t = true
	if (t) {
		call reuse_multiply as conditional_single {
			input: first = 1, second = 2
		}
	}

	# Frag applet via conditional, multiple call
	if (!t) {
		call reuse_multiply as conditional_multi_1 {
			input: first = 3, second = 4
		}
		call reuse_multiply as conditional_multi_2 {
			input: first = 5, second = 6
		}
	}

	# Frag applet via scatter, single call
	Array[String] scatter_single_arr = ["Hello", "World"]
	scatter (i in scatter_single_arr) {
		call reuse_print as scatter_single {
			input: str = i
		}
	}

	# Frag applet via scatter, multiple call
	Array[String] scatter_multi_arr = ["Hello", "Good day", "Ahoj", "Dobry den"]
	String punctuation = "!"
	scatter (i in scatter_multi_arr) {
		call reuse_print as scatter_multi_1 {
			input: str = i
		}
		call reuse_print as scatter_multi_2 {
			input: str = punctuation
		}
	}

	# Do minor edit unrelated to frag applets
	call reuse_print as ZZZ {
		input: str = "Unrelated"
	}
}