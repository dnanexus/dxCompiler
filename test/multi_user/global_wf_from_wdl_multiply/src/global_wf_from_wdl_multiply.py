#!/usr/bin/env python3

import dxpy

@dxpy.entry_point('main')
def main(multiply_first, multiply_second):
    product = multiply_first * multiply_second
    output = {}
    output["product"] = product
    return output

dxpy.run()
