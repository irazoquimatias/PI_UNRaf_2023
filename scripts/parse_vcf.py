#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import argparse

def parse_options():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", type=str, help="Input VCF file")
    parser.add_argument("-p", "--passed", type=str, help="Output VCF file with passed variants")
    parser.add_argument("-f", "--failed", type=str, help="Output VCF file with failed variants")
    parser.add_argument("-t", "--threshold", type=int, default = 20, help="Coverage threshold to discard variants")
    args = parser.parse_args()
    return args

def parse_vcf(args):
    fail_report = ""
    pass_report = ""
    with open(args.input) as i_vcf:
        for line in i_vcf:
            if line[0] == "#":
                fail_report += line
                pass_report += line
            else: 
                fields = line.split("\t")
                # Discard if there are frameshifts
                bases = len(fields[4]) - len(fields[3])
                if bases % 3 != 0:
                    fail_report += line
                else:
                    pass_report += line

                # Discard if the variant has low coverage
                info = fields[7].split(";")
                dp = int(info[0].split("=")[1])
                if dp < args.threshold:
                    fail_report += line
                else:
                    pass_report += line
                
                # Discard if low confidence
                sample = fields[9].split(":")
                gt = sample[0]
                if gt != "1/1":
                    fail_report += line
                else:
                    pass_report += line

    with open(args.passed, "w") as p_vcf, open(args.failed, "w") as f_vcf:
        p_vcf.write(pass_report)
        f_vcf.write(fail_report)

if __name__ == "__main__":
    args = parse_options()
    parse_vcf(args)

