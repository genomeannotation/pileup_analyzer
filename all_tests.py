#!/usr/bin/env python

import unittest
from pileup_parser_classes import ConsensusCaller, PileSanitizer

class TestPileSanitizer(unittest.TestCase):
    def test_sanitize(self):
        sanitizer = PileSanitizer()
        pile1 = 'AAAA+3CCCGG'
        result1 = sanitizer.sanitize(pile1)
        self.assertEqual('AAAAGG', result1)

class TestConsensusCaller(unittest.TestCase):
    def test_initialize(self):
        caller = ConsensusCaller(depth= 10, frequency = .99)
        self.assertEqual(10, caller.min_depth_of_coverage)

    def test_call(self):
        caller = ConsensusCaller(depth = 5, frequency = .9)
        test_pile1 = 'AAAAAAAAAAAAAAACA'
        test_pile2 = 'AAAAAAAAAAAAAAAAAAA'
        pile = [test_pile1, test_pile2]
        result = caller.call(pile)
        self.assertEqual('A', result)

##################################
if __name__ == '__main__':
    unittest.main()
