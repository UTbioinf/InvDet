#!/usr/bin/env python

import sys
import argparse
import os
import pysam

def run(args):
    fin_bam = pysam.AlignmentFile(args.bam, "rb")
    print os.path.join(args.directory + "simple_delta")
    fout = open(os.path.join(args.directory, "simple_delta"), "w") if args.directory else sys.stdout
    # write total # of references
    fout.write("{}\n".format( fin_bam.nreferences ))
    # write reference name: <length> <reference name>
    lengths = fin_bam.lengths
    for i in xrange(fin_bam.nreferences):
        fout.write("{} {}\n".format( lengths[i], fin_bam.get_reference_name(i) ))
    # write each alignment header
    if args.directory:
        fout.close()



def parse_args( argv = None ):
    parser = argparse.ArgumentParser(description = "Bam file extractor")
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    #parser.add_argument("-t", "--target", required=True, help="contig/scaffold file (fasta/q)")
    #parser.add_argument("-r", "--read", required=True, help="long reads file (fasta/q)")
    parser.add_argument("-d", "--directory", help="directory of the output file. If not set, it will be output to stdout")
    return parser.parse_args( argv )

def main( argv = None ):
    args = parse_args( argv )
    run( args )
    

if __name__ == "__main__":
    main(["-b", "/Volumes/bioinfo/tmp/finding_inversions/afun_out.ctg.bam", "-d", "."])
