#!/usr/bin/env python

import sys
import csv
from pileup_parser_classes import *

if len(sys.argv) != 2:
    sys.stderr.write("usage: python pileup_parser.py <input.pileup>\n")
    sys.exit()

tsv_file = sys.argv[1]
output_file = sys.argv[1]+".results"

## Create parser that identifies first 6 piles as control,
## next 6 as experimental 
## TODO take this info as cmdline arg?
groups = [[0,1,2,3,4,5], [6,7,8,9,10,11]]
parser = PileupLineParser(groups)

outfile = open(output_file, 'w')
writer = csv.writer(outfile, delimiter='\t', quoting=csv.QUOTE_NONE)

def write_error(line, error):   
    msg = "Discarded "+line[0]+", locus "+line[1]"; Reason: "+error
    sys.stderr.write(msg)

with open(tsv_file, 'rb') as file:
    reader = csv.reader(file, delimiter='\t', quotechar='|')
    for line in reader:
        ## validate with parser
        ## if not, write info to stderr, continue
        ## write_error(line, "Nominal depth below threshold")

        ## generate Locus with parser

        ## locus.sanitize_all
        ## locus.filter_all
        ## locus.validate_depth
        ## if not, write info to stderr, continue
        ## write_error(line, "Depth below threshold after quality filtering")

        ## locus.call_consensus
        ## if not, write info to stderr, continue
        ## write_error(line, "Unable to call consensus")

        ## locus.generate_stats(call)
        ## writer.writerow(stats)

outfile.close()

