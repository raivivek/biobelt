#! /usr/bin/env python3

import sys
import os
import argparse

import pysam

__description__="""
    Fetch cell barcode from 'CB' tag and append to read name
    for processing with snapATAC.
"""

# SnapATAC requires:
#  (a) barcode to be "prepended" to the read name, and
#  (b) name sorted

parser = argparse.ArgumentParser(description=__description__)
parser.add_argument('bam_in', help='Bam file to adjust read names in')
parser.add_argument('-o', help='Output file [default: STDOUT]',
                    required=False, default=sys.stdout)
args = parser.parse_args()

f = pysam.AlignmentFile(args.bam_in, 'rb')
f_new = pysam.AlignmentFile(args.o, 'wb', template=f)

count = 0
for read in f.fetch(until_eof=True):
    count += 1
    if count % 1000000 == 0:
        print(f'Processed {count} reads so far', file=sys.stderr)

    # parse out the barcode
    barcode = read.get_tag('CB')
    # update the barcode
    read.query_name = f"{barcode}:{read.query_name}"
    f_new.write(read)
