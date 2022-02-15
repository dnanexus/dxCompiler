version 1.1

struct TestStruct {
	File aa
	File? bb
	File cc
	File? dd
}

task apps1052_optional_block_inputs_wdl11_t1 {
	command <<<
		echo a >a
		echo b >b
		echo c >c
	>>>
	output {
		TestStruct testStructOut = TestStruct { 
			aa: "a",
			bb: "b",
			cc: "c",
			dd: "d"
		}
	}
}

workflow apps1052_optional_block_inputs_wdl11 {
	input {
	}
	call apps1052_optional_block_inputs_wdl11_t1 as t1
	File aa = t1.testStructOut.aa # T to T
	File? bb = t1.testStructOut.bb # T? to T?
	File? cc = t1.testStructOut.cc # T to T?
	File? dd = t1.testStructOut.dd # null to T?
	output {
		File aa_out = aa
		File? bb_out = bb
		File? cc_out = cc
		File? dd_out = dd
	}
}
