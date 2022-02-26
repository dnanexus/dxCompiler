version 1.0

struct innerStruct {
	Int? aa
	Int? bb
    }

struct TestStruct {
    innerStruct inner
}

workflow apps_1052_optional_compound_input_wdl10 {
    input {
        TestStruct t1
    }
    
    TestStruct t2 = {"inner": {"aa": 2, "bb": 3}}
    Int? a = select_first([t1.inner.aa, t1.inner.bb, t2.inner.aa, t2.inner.bb])

    output {
        Int? out_int = a
    }
}
