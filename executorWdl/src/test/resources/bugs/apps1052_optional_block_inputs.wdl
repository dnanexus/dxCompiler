version 1.1

struct TestStruct {
	File aa
	File? bb
}

task apps1052_optional_block_inputs_t1 {
	command <<<
		echo a >a
		echo b >b
	>>>
	output {
		TestStruct testStructOut = TestStruct { 
			aa: "a",
			bb: "b"
		}
	}
}

workflow apps1052_optional_block_inputs {
	input {
	}
	call apps1052_optional_block_inputs_t1 as t1
	File aa = t1.testStructOut.aa # T to T
	File? bb = t1.testStructOut.bb # T? to T?
	output {
		File aa_out = aa
		File? bb_out = bb
	}
}
