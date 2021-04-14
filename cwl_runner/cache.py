#!/usr/bin/env python3


cache = {}


def get_cache(location):
    if location in cache:
        return cache[location]


def add_to_cache(local, remote):
    cache[local] = remote
