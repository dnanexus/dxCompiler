def rm_suffix(original_string: str, suffix: str) -> str:
    if not original_string:
        return ""
    if not suffix:
        return original_string
    if original_string[-1] == suffix[-1]:
        return rm_suffix(original_string[:-1], suffix[:-1])
    else:
        raise ValueError("rm_suffix(): suffix is not present in the original string")


def rm_prefix(original_string: str, prefix: str) -> str:
    if not original_string:
        return ""
    if not prefix:
        return original_string
    if original_string[0] == prefix[0]:
        return rm_prefix(original_string[1:], prefix[1:])
    else:
        raise ValueError("rm_suffix(): prefix is not present in the original string")