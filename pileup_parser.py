#!/usr/bin/env python

import sys
import csv

if len(sys.argv) != 2:
    sys.stderr.write("usage: python do_county_thing_with_pileup.py <input.pileup>\n")
    sys.exit()

tsv_file = sys.argv[1]

with open(tsv_file, 'rb') as file:
    reader = csv.reader(file, delimiter='\t', quotechar='|')
    for line in reader:
        print(str(line))

"""sample_output = []
sample_output.append('writer writes lists...')
sample_output.append('like this one.')

writer = csv.writer(sys.stdout, delimiter='\t', quoting=csv.QUOTE_NONE)
writer.writerow(sample_output)"""
