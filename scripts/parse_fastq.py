#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import matplotlib.pyplot as plt
import argparse
import os
import re

def parse_options():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", type=str, help="Input folder with unfiltered FASTQ files")
    parser.add_argument("-o", "--output", type=str, help="Output folder for filtered FASTQ files")
    parser.add_argument("--min_q", type=int, default = 12, help="Min quality to accept a read (default: 12)")
    parser.add_argument("--min_len", type=int, default = 1200, help="Min length to accept a read (default: 1200)")
    parser.add_argument("--max_len", type=int, default = 1700, help="Max length to accept a read (default: 1700)")
    parser.add_argument("-p", "--plot", type=str, help="Plot bar chart of accepted/rejected reads (default: don't do it)")
    parser.add_argument("-v", '--verbose', action="store_true")
    args = parser.parse_args()
    return args

def plot_reads(args,files, rejected, accepted):
    plt.bar(files, accepted, color='#428bca')
    plt.bar(files, rejected, bottom=accepted, color='#d9534f')
    plt.legend(["Retained", "Discarded"], loc =(1.04, 0.5))
    plt.xticks(rotation = 90)
    plt.tight_layout()
    plt.savefig(args.plot)

def parse_fastq (args):
    map_reads_barcode = {}
    files = []
    rejected = []
    accepted = []
    for filename in os.listdir(args.input):
        print ("Working on file", filename)  if args.verbose else None
        fullname = os.path.join(args.input, filename)
        line_number = 0
        read = ""
        report = ""
        accepted_reads = 0
        total_reads = 0
        id_read = ""
        with open(fullname, "r") as ifh:
            for line in ifh:
                line_number += 1
                read += line
                if line_number%4 == 0:
                    total_reads += 1
                    if len(line) > args.min_len and len(line) < args.max_len:
                        line = line.rstrip()
                        quality = 0
                        for c in line:
                            quality += (ord(c) - 33)
                        if quality/len(line) > args.min_q:
                            report += read
                            accepted_reads += 1
                            map_reads_barcode[id_read]={"barcode":filename}
                    read = ""
                elif (line_number-1) %4 == 0:
                    id_read = line.split(" ")[0]
                    id_read = id_read[1:]

        fullname = os.path.join(args.output, filename)
        with open(fullname, "w") as ofh:
            ofh.write(report)

        files.append(filename)
        rejected.append(total_reads-accepted_reads)
        accepted.append(accepted_reads)

    if args.plot:
        print ("Ploting") if args.verbose else None
        plot_reads(args,files, rejected, accepted)

    return map_reads_barcode

if __name__ == "__main__":
    args = parse_options()
    if not re.match('^/', args.output):
        args.output =  os.getcwd() + '/' + args.output
    if os.path.exists(args.output) and os.listdir(args.output):
        print('ERROR: la carpeta de salida no esta vacia')
        exit()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    read_barcode = parse_fastq(args)

