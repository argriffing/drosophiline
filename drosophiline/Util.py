
import unittest

import iterutils

def stripped_lines(lines):
    """
    This function yields nonempty stripped lines.
    """
    for line in lines:
        line = line.strip()
        if line:
            yield line

def get_stripped_lines(lines):
    """
    @return: a list of nonempty stripped lines
    """
    return list(iterutils.stripped_lines(lines))

def hamming_distance(first, second):
    return sum(1 for a, b in zip(first, second) if a != b)

class Cache:
    def __init__(self, callback, cache_limit):
        self.callback = callback
        self.cache_limit = cache_limit
        self.cache = {}
    def __call__(self, arg):
        try:
            return self.cache[arg]
        except KeyError:
            value = self.callback(arg)
            if not self.is_full():
                self.cache[arg] = value
            return value
    def is_full(self):
        if self.cache_limit is None:
            return False
        return (len(self.cache) >= self.cache_limit)


class TestUtil(unittest.TestCase):

    def test_placeholder(self):
        pass


if __name__ == '__main__':
    unittest.main()
