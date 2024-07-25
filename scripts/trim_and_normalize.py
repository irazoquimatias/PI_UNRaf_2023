#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import pysam
import argparse

def parse_options():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", type=str, help="Input folder with unfiltered FASTQ files")
    parser.add_argument("-o", "--output", type=str, help="Output folder for filtered FASTQ files")
    parser.add_argument("-p","--primers", type=str, help="BED file with primers positions")
    parser.add_argument("-v", '--verbose', action="store_true")
    args = parser.parse_args()
    return args

def pos_in_primer(primers, pos):
    output = None
    for p in primers:
        if pos >= (primers[p]["start"]-1) and pos < (primers[p]["end"]-1):
            output = primers[p]["start"]-1,primers[p]["end"]-1
    return output

def load_bed_file(args):
    primers = {}
    with open (args.primers) as ifh:
        for l in ifh:
            l = l.rstrip()
            fields = l.split("\t")
            if int(fields[1]) < int(fields[2]):
                primers[fields[3]] = {"start": int(fields[1]), "end": int(fields[2])}
            else:
                primers[fields[3]] = {"start": int(fields[2]), "end": int(fields[1])}
    return primers

def length_cigar(cigar):
    length = 0
    for t in cigar:
        if t[0] == 0 or t[0] == 1 or t[0] == 4:
            length += t[1]
    return length

def soft_clip_primers(args, primers):
    check = "NADA"

    for read in bf_in.fetch():
        f_offset = 0
        r_offset = 0
        if not read.is_unmapped:
            pos = pos_in_primer(primers, read.reference_start)
            if pos:
                f_offset = pos[1] - read.reference_start + 1
            pos = pos_in_primer(primers, read.reference_end)
            if pos:
                r_offset = read.reference_end - pos[0]
        
        if read.query_name == check:
            print(f_offset, r_offset, read.query_length)
            print("pre", read.cigar)

        if f_offset > 0:
            offset_tuple=(4, f_offset)
            first_tuple=(None, None)
            extra_tuple=(None, None)
            while f_offset > 0:
                cigar_tuple = read.cigar[0]
                read.cigar = read.cigar[1:]

                print ("IN_pre", f_offset, cigar_tuple) if read.query_name == check else 0
                if cigar_tuple[0] == 1:   # If there's an insertion in the primer region, increase the sequence to clip
                    offset_tuple = (offset_tuple[0], offset_tuple[1] + cigar_tuple[1])
                elif cigar_tuple[0] == 5 and cigar_tuple[1] >= f_offset: # If there's already a clip in the region there's no need to clip further
                    offset_tuple = (5, offset_tuple[1] + cigar_tuple[1] - f_offset)
                    f_offset = 0
                elif cigar_tuple[0] == 5 and cigar_tuple[1] < f_offset: # If there's a region smaller than the primer region but a hard clip
                    first_tuple = cigar_tuple
                    f_offset -= cigar_tuple[1] 
                    offset_tuple = (4, f_offset)
                elif cigar_tuple[0] == 4 and cigar_tuple[1] >= f_offset: # If there's already a clip in the region there's no need to clip further
                    offset_tuple = (4, offset_tuple[1] + cigar_tuple[1] - f_offset)
                    f_offset = 0 
                elif cigar_tuple[0] == 4 and cigar_tuple[1] < f_offset: # If there's already a clip in the region there's no need to clip further
                    f_offset -= cigar_tuple[1] 
                elif cigar_tuple[0] == 0 and cigar_tuple[1] >= f_offset: # If there's a mapped region
                    extra_tuple = (0, cigar_tuple[1]-f_offset)
                    f_offset = 0 
                elif cigar_tuple[0] == 0 and cigar_tuple[1] < f_offset: # If there's a region smaller than the primer region
                    f_offset -= cigar_tuple[1] 
                print ("IN_post", f_offset, cigar_tuple) if read.query_name == check else 0

            if extra_tuple[0] == 0 and first_tuple[0] == 5:
                read.cigar = [first_tuple, offset_tuple, extra_tuple] + read.cigar
            elif extra_tuple[0] == 0:
                read.cigar = [offset_tuple, extra_tuple] + read.cigar
            elif first_tuple[0] == 5:
                read.cigar = [first_tuple, offset_tuple] + read.cigar
            else:
                read.cigar = [offset_tuple] + read.cigar

        if r_offset > 0:
            offset_tuple=(4, r_offset)
            extra_tuple=(None, None)
            last_tuple=(None, None)
            while r_offset > 0:
                cigar_tuple = read.cigar[-1]
                read.cigar = read.cigar[:-1]

                print ("IN_pre", r_offset, cigar_tuple) if read.query_name == check else 0
                if cigar_tuple[0] == 1:   # If there's an insertion in the primer region, increase the sequence to clip
                    offset_tuple = (offset_tuple[0], offset_tuple[1] + cigar_tuple[1])
                elif cigar_tuple[0] == 5 and cigar_tuple[1] >= r_offset: # If there's already a clip in the region there's no need to clip further
                    offset_tuple = (5, offset_tuple[1] + cigar_tuple[1] - r_offset)
                    r_offset = 0
                elif cigar_tuple[0] == 5 and cigar_tuple[1] < r_offset: # If there's a region smaller than the primer region but a hard clip
                    last_tuple = cigar_tuple
                    r_offset -= cigar_tuple[1] 
                    offset_tuple = (4, r_offset)
                elif cigar_tuple[0] == 4 and cigar_tuple[1] >= r_offset: # If there's already a clip in the region there's no need to clip further
                    offset_tuple = (4, offset_tuple[1] + cigar_tuple[1] - r_offset)
                    r_offset = 0 
                elif cigar_tuple[0] == 4 and cigar_tuple[1] < r_offset: # If there's a region smaller than the primer region but a hard clip
                    r_offset -= cigar_tuple[1] 
                elif cigar_tuple[0] == 0 and cigar_tuple[1] >= r_offset: # If there's a mapped region, clip the primer and modify the existing mapped region
                    extra_tuple = (0, cigar_tuple[1]-r_offset)
                    r_offset = 0 
                elif cigar_tuple[0] == 0 and cigar_tuple[1] < r_offset: # If there's a region smaller than the primer region, clip the whole region
                    r_offset -= cigar_tuple[1] 
                print ("IN_post", r_offset, cigar_tuple) if read.query_name == check else 0

            if extra_tuple[0] == 0 and last_tuple[0] == 5:
                read.cigar = read.cigar + [extra_tuple, offset_tuple, last_tuple]
            elif extra_tuple[0] == 0:
                read.cigar = read.cigar + [extra_tuple, offset_tuple]
            elif extra_tuple[0] == 0:
                read.cigar = read.cigar + [offset_tuple, last_tuple]
            else:
                read.cigar = read.cigar + [offset_tuple]

        if read.query_name == check:
            print("post", read.cigar)

        bf_out.write(read)

if __name__ == "__main__":
    args = parse_options()
    primers = load_bed_file(args)
    bf_in = pysam.AlignmentFile(args.input, "rb")
    bf_out = pysam.AlignmentFile(args.output, "wb", template=bf_in)
    soft_clip_primers(args, primers)
    bf_in.close()
    bf_out.close()

