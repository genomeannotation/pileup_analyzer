#!/usr/bin/env python

import unittest
from pileup_parser_classes import ConsensusCaller, PileSanitizer, QualityFilter

class TestPileSanitizer(unittest.TestCase):
    def test_sanitize(self):
        sanitizer = PileSanitizer()
        # remove insertions
        pile1 = 'AAAA+3CCCGG'
        result1 = sanitizer.sanitize(pile1)
        self.assertEqual('AAAAGG', result1)
        # remove deletions
        pile2 = 'AA-2CCAAGG'
        result2 = sanitizer.sanitize(pile2)
        self.assertEqual('AAAAGG', result2)
        # remove those '^' CIGAR thingies
        pile3 = 'AA^#AAGG'
        result3 = sanitizer.sanitize(pile3)
        self.assertEqual('AAAAGG', result3)
        # remove those '$' CIGAR thingies
        pile4 = 'A$(AAAGG'
        result4 = sanitizer.sanitize(pile4)
        self.assertEqual('AAAAGG', result4)

class TestQualityFilter(unittest.TestCase):
    def test_initialize(self):
        fil = QualityFilter(30)
        self.assertEqual(30, fil.quality_threshold)

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
