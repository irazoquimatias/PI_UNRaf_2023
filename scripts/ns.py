#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import Bio.SeqIO as bpio
import sys

for record in bpio.parse(sys.argv[1], "fasta"):
    seq = str(record.seq)
    length = len(record.seq)
    length_no_n = len(seq.strip("N"))
    ns = seq.count("N")
    ns_perc = ns/length*100
    ns_perc = "{:.2f}".format(ns_perc)
    print (record.id,"\t",length,"\t",length_no_n,"\t",ns,"\t",ns_perc)

