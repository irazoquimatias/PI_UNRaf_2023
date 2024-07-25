#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import pysam
import itertools
import argparse

def parse_options():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", type=str, help="Input BAM against reference")
    parser.add_argument("-o", "--output", type=str, help="Output BED file")
    parser.add_argument("-c","--coverage", type=int, default = 20, help="Min coverage to not be masked")
    parser.add_argument("-v", '--verbose', action="store_true")
    #TODO: add verbose
    args = parser.parse_args()
    return args

# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def eval_depth(args):
    # Array to save all positions with depth lower than threshold
    to_mask = []
    # Read file in bam format
    bf = pysam.AlignmentFile(args.input, "rb")
    reference = bf.references[0]
    # Sum all reads in each column in reference
    for pileupcolumn in bf.pileup(bf.references[0]):
        depth = 0
        for read in pileupcolumn.pileups:
            depth+=1
        to_mask.append(pileupcolumn.pos + 1) if depth < args.coverage else None
    # Convert list of numbers to intervals
    intervals = list(intervals_extract(to_mask))
    bf.close()
    return intervals, reference

def print_outfile(args, intervals, reference):
    report = ""
    for i in intervals:
        report = report + reference + "\t" + str(i[0])  + "\t" + str(i[1]) + "\n"
    with open (args.output, "w") as OFH:
        OFH.write(report)
    
if __name__ == "__main__":
    args = parse_options()
    intervals, reference = eval_depth(args)
    print_outfile(args, intervals, reference)

