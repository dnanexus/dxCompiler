import pytest
from dxcint.utils import rm_suffix, rm_prefix


def test_rm_suffix():
    assert rm_suffix("foo.cwl", ".cwl") == "foo"
    assert rm_suffix("foo.cwl", "") == "foo.cwl"
    assert rm_suffix("", "foo") == ""
    with pytest.raises(ValueError) as e:
        assert rm_suffix("foo.cwl", "foo")
        assert e


def test_rm_prefix():
    assert rm_prefix("foo.cwl", "foo") == ".cwl"
    assert rm_prefix("foo.cwl", "") == "foo.cwl"
    assert rm_prefix("", "foo") == ""
    with pytest.raises(ValueError) as e:
        assert rm_prefix("bar.cwl", "foo")
        assert e
