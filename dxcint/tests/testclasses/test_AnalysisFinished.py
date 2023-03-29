def test__validate(build_executable_wdl):
    result = build_executable_wdl.get_test_result()
    assert result.get("passed")
