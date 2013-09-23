#!/usr/bin/env python

import re

class PileSanitizer:
    def __init__(self):
        self.indel_re = re.compile('[\+\-]')
        self.cigar_re = re.compile('[\^\$]')

    def sanitize(self, pile):
        clean_pile = ''
        skip = 0
        on_indel = False
        for base in pile:
            m1 = self.indel_re.match(base)
            if m1:
                on_indel = True
                continue
            m2 = self.cigar_re.match(base)
            if m2:
                skip = 1
                continue
            if on_indel:
                skip = int(base)
                on_indel = False
                continue
            if skip > 0:
                skip -= 1
                continue
            clean_pile += base
        return clean_pile

class QualityFilter:
    def __init__(self, min_qual):
        self.quality_threshold = min_qual

class ConsensusCaller:

    def __init__(self, depth, frequency):
        self.min_depth_of_coverage = depth
        self.min_base_frequency = frequency

    # calls consensus on list of piles
    def call(self, piles):
        total_length = 0
        counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
        for pile in piles:
            if len(pile) < self.min_depth_of_coverage:
                return None
            total_length += len(pile)
            for base in pile:
                if base == 'A' or base == 'a':
                    counts['A'] += 1
                elif base == 'C' or base == 'c':
                    counts['C'] += 1
                elif base == 'T' or base == 't':
                    counts['T'] += 1
                elif base == 'G' or base == 'g':
                    counts['G'] += 1
        ref_base = max(counts.iterkeys(), key=(lambda key: counts[key]))
        if (float(counts[ref_base]) / total_length) >= self.min_base_frequency:
            return ref_base
        else:
            return None
