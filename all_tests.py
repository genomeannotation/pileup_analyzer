#!/usr/bin/env python

import unittest
from pileup_parser_classes import ConsensusCaller

class TestConsensusCaller(unittest.TestCase):
    def test_initialize(self):
        caller = ConsensusCaller(depth= 10, frequency = .99)
        self.assertEqual(10, caller.min_depth_of_coverage)



##################################
if __name__ == '__main__':
    unittest.main()
