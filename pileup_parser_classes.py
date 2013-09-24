#!/usr/bin/env python

import re
from phred import PhredHelper

class Pile:
    def __init__(self, bases, scores):
        self.bases = bases
        self.scores = scores

# Operates on a Pile.bases string
class PileSanitizer:
    def __init__(self):
        self.indel_re = re.compile('[\+\-]')
        self.cigar_re = re.compile('[\^\$]')

    def sanitize(self, bases):
        clean_bases = ''
        skip = 0
        on_indel = False
        for base in bases:
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
            clean_bases += base
        return clean_bases

# Operates on a Pile object
class QualityFilter:
    def __init__(self, min_qual, offset=33):
        self.quality_threshold = min_qual
        self.phred_helper = PhredHelper(offset)

    def filter(self, pile):
        keep_indices = []
        keep_bases = ''
        keep_scores = ''
        # build list of indices of bases/scores to keep
        for index, score in enumerate(pile.scores):
            if self.phred_helper.char_to_int(score) >= self.quality_threshold:
                keep_indices.append(index)
        # build new bases and scores strings
        for n in keep_indices:
            keep_bases += pile.bases[n]
            keep_scores += pile.scores[n]
        # update pile fields
        pile.bases = keep_bases
        pile.scores = keep_scores

# Operates on a list of strings (each string is a Pile.bases)
class ConsensusCaller:
    def __init__(self, depth, frequency):
        self.min_depth_of_coverage = depth
        self.min_base_frequency = frequency

    # calls consensus on list of pile.bases strings
    def call(self, list_of_bases):
        total_length = 0
        counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
        for pile in list_of_bases:
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

# Takes a list of two lists (0-indexed sample numbers for control and experimental groups)
class PileupLineParser:
    def __init__(self, grouplist):
        if len(grouplist) != 2:
            sys.err.write("PileupLineParser: grouplist must be a list of 2 lists.")
        self.control_group = grouplist[0]
        self.experimental_group = grouplist[1]

    def get_piles_from_group(self, group, line):
        piles = []
        for i in group:
            starting_index = 3 * (i+1)
            bases = line[starting_index + 1]
            scores = line[starting_index + 2]
            piles.append(Pile(bases, scores))
        return piles

    def get_control_piles(self, line):
        return self.get_piles_from_group(self.control_group, line)

    def get_experimental_piles(self, line):
        return self.get_piles_from_group(self.experimental_group, line)

    def get_all_bases(self, line):
        all_bases = []
        number_of_samples = (len(line) - 3) / 3
        for i in range(number_of_samples):
            starting_index = 3 * (i+1)
            bases = line[starting_index + 1]
            all_bases.append(bases)
        return all_bases

    def get_lengths(self, line):
        lengths = []
        number_of_samples = (len(line) - 3) / 3
        for i in range(number_of_samples):
            length_index = 3 * (i+1)
            lengths.append(int(line[length_index]))
        return lengths

    def validate(self, line, min_length):
        lengths = self.get_lengths(line)
        for length in lengths:
            if length < min_length:
                return False
        return True

    def get_chromosome(self, line):
        return line[0]

    def get_coordinate(self, line):
        return line[1]

    def generate_locus(self, line):
        chrom = self.get_chromosome(line)
        coord = self.get_coordinate(line)
        ctrl = self.get_control_piles(line)
        exp = self.get_experimental_piles(line)
        return Locus(chrom, coord, ctrl, exp)

class Locus:
    def __init__(self, chrom, coord, ctrl_piles, exp_piles):
        self.chromosome = chrom
        self.coordinate = coord
        self.control_piles = ctrl_piles
        self.experimental_piles = exp_piles

    def to_string(self):
        result = "chromosome: " + self.chromosome + "; "
        result += "coordinate: " + self.coordinate + "; "
        result += "control piles: " + str(len(self.control_piles)) + "; "
        result += "experimental piles: " + str(len(self.experimental_piles))
        return result

    def sanitize_all(self):
        sani = PileSanitizer()        
        for pile in self.control_piles:
            pile.bases = sani.sanitize(pile.bases)
        for pile in self.experimental_piles:
            pile.bases = sani.sanitize(pile.bases)

    def filter_all(self, qscore, offset):
        filt= QualityFilter(qscore, offset)
        for pile in self.control_piles:
            filt.filter(pile)
        for pile in self.experimental_piles:
            filt.filter(pile)

    def validate_depth(self, depth):
        for pile in self.control_piles:
            if len(pile.scores) < depth:
                return False
        for pile in self.experimental_piles:
            if len(pile.scores) < depth:
                return False
        return True
        
        




