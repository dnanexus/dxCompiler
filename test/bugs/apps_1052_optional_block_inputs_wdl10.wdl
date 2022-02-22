version 1.0

struct TestStruct {
	File aa
	File? bb
	File cc
	File? dd
}

task apps_1052_optional_block_inputs_wdl10_t1 {
	command <<<
		echo a >a
		echo b >b
		echo c >c
	>>>
	output {
		TestStruct testStructOut = object { 
			aa: "a",
			bb: "b",
			cc: "c",
			dd: "d"
		}
	}
}

workflow apps_1052_optional_block_inputs_wdl10 {
	input {
	}
	call apps_1052_optional_block_inputs_wdl10_t1 as t1
	File aa = t1.testStructOut.aa # T to T
	File? bb = t1.testStructOut.bb # T? to T?
	File? cc = t1.testStructOut.cc # T to T?
	File? dd = t1.testStructOut.dd # null to T?
	output {
		# Expect 3 files to exist
		Int filesCount = length(select_all([aa, bb, cc, dd]))
	}
}
