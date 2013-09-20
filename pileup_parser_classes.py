#!/usr/bin/env python

class ConsensusCaller:

    def __init__(self, depth, frequency):
        self.min_depth_of_coverage = depth
        self.min_base_frequency = frequency
