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
        # remove those '^' CIGAR thingies and the base after
        bases3 = 'AA^#AAGG'
        result3 = sanitizer.sanitize(bases3)
        self.assertEqual('AAAAGG', result3)
        # remove those '$' CIGAR thingies
        bases4 = 'A$AAAGG'
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
        caller = ConsensusCaller(frequency = .99)
        self.assertEqual(.99, caller.min_base_frequency)

    def test_call(self):
        caller = ConsensusCaller(frequency = .9)
        test_bases1 = 'AAAAAAAAAAAAAAACA'
        test_bases2 = 'AAAAAAAAAAAAAAAAAAA'
        bases = [test_bases1, test_bases2]
        result = caller.call(bases)
        self.assertEqual('A', result)

class TestPileupLineParser(unittest.TestCase):
    def setUp(self):
        groups = [[0, 1], [2, 3]]
        self.parser = PileupLineParser(groups)
        input_string = """comp102583_c0_seq1      667     N       34      AAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa      B::D0DDBBDBDDDDDD3DDD5DDD>@DBDDDDD      14      AAAAAAAaaaaaaa  DDDDDBBDDDDDDD  40      AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaa        DD65DDBDDDDDDDDD86BB#DDDDDDDBDDDDD@DD;D9        30      AAAAAAAAAAAAAAAaaaaaaacaaaaaaa  DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD"""
        self.test_input = input_string.split()
        
    def test_initialize(self):
        self.assertEqual([0, 1], self.parser.control_group)

    def test_get_control_piles(self):
        control_piles = self.parser.get_control_piles(self.test_input)
        self.assertEqual(2, len(control_piles))
        self.assertEqual('AAAAAAAaaaaaaa', control_piles[1].bases)

    def test_get_experimental_piles(self):
        exp_piles = self.parser.get_experimental_piles(self.test_input)
        self.assertEqual(2, len(exp_piles))
        self.assertEqual('DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD', exp_piles[1].scores)

    def test_get_all_bases(self):
        all_bases = self.parser.get_all_bases(self.test_input)
        self.assertEqual(4, len(all_bases))
        self.assertEqual('AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaa', all_bases[2])

    def test_get_lengths(self):
        lengths = self.parser.get_lengths(self.test_input)
        expected = [34, 14, 40, 30]
        self.assertEqual(expected, lengths)

    def test_validate(self):
        self.assertTrue(self.parser.validate(self.test_input, 10))
        self.assertFalse(self.parser.validate(self.test_input, 15))
        
    def test_get_chromosome(self):
        self.assertEqual('comp102583_c0_seq1', self.parser.get_chromosome(self.test_input))

    def test_get_coordinate(self):
        self.assertEqual('667', self.parser.get_coordinate(self.test_input))

    def test_generate_locus(self):
        locus = self.parser.generate_locus(self.test_input)
        self.assertEqual('Locus', locus.__class__.__name__)   #assertIsInstance only in py2.7!
        self.assertEqual('DDDDDBBDDDDDDD', locus.control_piles[1].scores)
        
class TestLocus(unittest.TestCase):
    def setUp(self):
        chrom = 'comp102583_c0_seq1'
        coord = '667'
        pile1 = Pile('AAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa', 'B::D0DDBBDBDDDDDD3DDD5DDD>@DBDDDDD')
        pile2 = Pile('AAAAAAAaaaaaaa', 'DDDDDBBDDDDDDD')
        control_piles = [pile1, pile2]
        pile3 = Pile('AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaa', 'DD65DDBDDDDDDDDD86BB#DDDDDDDBDDDDD@DD;D9')
        pile4 = Pile('AAAAAAAAAAAAAAAaaaaaaacaaaaaaa', 'DD6DD@DDBDD8DBDDDD5DDD#DDDBDDD')
        exp_piles = [pile3, pile4]
        self.locus = Locus(chrom, coord, control_piles, exp_piles)

    def setUp2(self):
        chrom = 'test_chrom'
        coord = '123'
        pile1 = Pile('GAT+3CCGTACA', 'hghhghg')
        pile2 = Pile('GATTA^&CA', 'hghhghg')
        control_piles = [pile1, pile2]
        pile3 = Pile('GATTAC-2GGA', 'hghhghg')
        pile4 = Pile('G$ATTA$CA', 'hghhghg')
        exp_piles = [pile3, pile4]
        self.locus2 = Locus(chrom, coord, control_piles, exp_piles)

    def setUp3(self):
        chrom = 'test_chrom'
        coord = '123'
        pile1 = Pile('GATTACA', 'hghhghg')
        pile2 = Pile('GATTACA', 'ggggggg')
        control_piles = [pile1, pile2]
        pile3 = Pile('GATTACA', 'hijklmn')
        pile4 = Pile('GATTACA', 'hghhghg')
        exp_piles = [pile3, pile4]
        self.locus3 = Locus(chrom, coord, control_piles, exp_piles)

    def test_init(self):
        self.assertEqual('comp102583_c0_seq1', self.locus.chromosome)

    def test_to_string(self):
        expected = "chromosome: comp102583_c0_seq1; coordinate: 667; "
        expected += "control piles: 2; experimental piles: 2"
        self.assertEqual(expected, self.locus.to_string())

    def test_sanitize_all(self):
        self.setUp2()
        self.locus2.sanitize_all()
        self.assertEqual('GATTACA', self.locus2.control_piles[0].bases)
        self.assertEqual('GATTACA', self.locus2.control_piles[1].bases)
        self.assertEqual('GATTACA', self.locus2.experimental_piles[0].bases)
        self.assertEqual('GATTACA', self.locus2.experimental_piles[1].bases)

    def test_filter_all(self):
        self.setUp3()
        self.locus3.filter_all(40, 64)
        self.assertEqual('GTTC', self.locus3.control_piles[0].bases)
        self.assertEqual('', self.locus3.control_piles[1].bases)
        self.assertEqual('GATTACA', self.locus3.experimental_piles[0].bases)

    def test_validate_depth(self):
        self.assertTrue(self.locus.validate_depth(14))
        self.assertFalse(self.locus.validate_depth(15))

    def test_call_consensus(self):
        call = self.locus.call_consensus(frequency = 0.9)
        self.assertEqual('A', call)
        self.setUp3()
        call2 = self.locus3.call_consensus(frequency = 0.42)
        self.assertEqual('A', call2)
        call3 = self.locus3.call_consensus(frequency = 0.43)
        self.assertFalse(call3)

    def test_generate_stats(self):
        expected = ['comp102583_c0_seq1_locus667']
        expected.extend([48, 48, 1.0, '|', 69, 70, 0.98571428571428577, '|'])
        expected.extend([34, 34, 14, 14, '|', 40, 40, 29, 30])
        self.assertEqual(expected, self.locus.generate_stats('A'))
        

        

        

        

##################################
if __name__ == '__main__':
    unittest.main()
