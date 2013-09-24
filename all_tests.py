#!/usr/bin/env python

import unittest
#from pileup_parser_classes import Pile, ConsensusCaller, PileSanitizer, QualityFilter, PileupLineParser
from pileup_parser_classes import *

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
        bases1 = 'AAAA+3CCCGG'
        result1 = sanitizer.sanitize(bases1)
        self.assertEqual('AAAAGG', result1)
        # remove deletions
        bases2 = 'AA-2CCAAGG'
        result2 = sanitizer.sanitize(bases2)
        self.assertEqual('AAAAGG', result2)
        # remove those '^' CIGAR thingies
        bases3 = 'AA^#AAGG'
        result3 = sanitizer.sanitize(bases3)
        self.assertEqual('AAAAGG', result3)
        # remove those '$' CIGAR thingies
        bases4 = 'A$(AAAGG'
        result4 = sanitizer.sanitize(bases4)
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
        test_bases1 = 'AAAAAAAAAAAAAAACA'
        test_bases2 = 'AAAAAAAAAAAAAAAAAAA'
        bases = [test_bases1, test_bases2]
        result = caller.call(bases)
        self.assertEqual('A', result)

class TestPileupLineParser(unittest.TestCase):
    def test_initialize(self):
        groups = [[0, 1], [2, 3]]
        parser = PileupLineParser(groups)
        self.assertEqual([0, 1], parser.control_group)

    def test_get_control_piles(self):
        groups = [[0, 1], [2, 3]]
        parser = PileupLineParser(groups)
        input_string = """comp102583_c0_seq1      667     N       34      AAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa      B::D0DDBBDBDDDDDD3DDD5DDD>@DBDDDDD      14      AAAAAAAaaaaaaa  DDDDDBBDDDDDDD  40      AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaa        DD65DDBDDDDDDDDD86BB#DDDDDDDBDDDDD@DD;D9        30      AAAAAAAAAAAAAAAaaaaaaacaaaaaaa  DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD"""
        test_input = input_string.split()
        control_piles = parser.get_control_piles(test_input)
        self.assertEqual(2, len(control_piles))
        self.assertEqual('AAAAAAAaaaaaaa', control_piles[1].bases)

    def test_get_experimental_piles(self):
        groups = [[0, 1], [2, 3]]
        parser = PileupLineParser(groups)
        input_string = """comp102583_c0_seq1      667     N       34      AAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa      B::D0DDBBDBDDDDDD3DDD5DDD>@DBDDDDD      14      AAAAAAAaaaaaaa  DDDDDBBDDDDDDD  40      AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaa        DD65DDBDDDDDDDDD86BB#DDDDDDDBDDDDD@DD;D9        30      AAAAAAAAAAAAAAAaaaaaaacaaaaaaa  DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD"""
        test_input = input_string.split()
        exp_piles = parser.get_experimental_piles(test_input)
        self.assertEqual(2, len(exp_piles))
        self.assertEqual('DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD', exp_piles[1].scores)
        
        

##################################
if __name__ == '__main__':
    unittest.main()
