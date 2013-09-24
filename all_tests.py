#!/usr/bin/env python

import unittest
from pileup_parser_classes import Pile, ConsensusCaller, PileSanitizer, QualityFilter

class TestPile(unittest.TestCase):
    def test_initialize(self):
        bases1 = 'GATTACA'
        scores1 = 'hhakasf'
        test_pile = Pile(bases1, scores1)
        self.assertEqual(bases1, test_pile.bases)
        self.assertEqual(scores1, test_pile.scores)

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
        fil = QualityFilter(30, 64)
        self.assertEqual(30, fil.quality_threshold)
        self.assertEqual(64, fil.phred_helper.offset)

    def test_filter(self):
        bases1 = 'GATTACA'
        scores1 = 'hUqpVu^'
        pile1 = Pile(bases1, scores1)
        fil = QualityFilter(40, 64)
        fil.filter(pile1)
        self.assertEqual('GTTC', pile1.bases)
        self.assertEqual('hqpu', pile1.scores)
        

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
