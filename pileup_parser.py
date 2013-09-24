#!/usr/bin/env python

import sys
import csv
from pileup_parser_classes import *

if len(sys.argv) != 2:
    sys.stderr.write("usage: python pileup_parser.py <input.pileup>\n")
    sys.exit()

tsv_file = sys.argv[1]
output_file = sys.argv[1]+".results"
groups = [[0,1,2,3,4,5], [6,7,8,9,10,11]]
min_depth = 1
min_qual = 1
offset = 33     # set to 33 for sanger, 64 for illumina 

## TODO take this info as cmdline arg?
parser = PileupLineParser(groups)

outfile = open(output_file, 'w')
writer = csv.writer(outfile, delimiter='\t', quoting=csv.QUOTE_NONE)

def write_error(line, error):   
    msg = "Discarded " + line[0] + ", locus " + line[1] + "; Reason: " + error
    sys.stderr.write(msg)

def write_header():
    # TODO should read groups so not hard-coded
    header = ["chromosome_locus", "control match", "control total", "|"]
    header.extend(["exp. match", "exp. total", "|"])
    header.extend(["each control sample match, total", "|", "each exp. sample match, total"])

with open(tsv_file, 'rb') as file:
    reader = csv.reader(file, delimiter='\t', quotechar='|')
    for line in reader:
        print(str(line))
        if not parser.validate(line, min_depth):
            write_error(line, "Nominal depth below threshold")
        else:
            locus = parser.generate_locus(line)
            locus.sanitize_all()
            locus.filter_all(min_qual, offset)
            if not locus.validate_depth(min_depth):
                write_error(line, "Depth below threshold after quality filtering")
            else:        
                call = locus.call_consensus
                if not call:
                    write_error(line, "Unable to call consensus")
                else:
                    stats = locus.generate_stats(call)
                    writer.writerow(stats)

outfile.close()

