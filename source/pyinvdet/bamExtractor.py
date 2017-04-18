#!/usr/bin/env python

import sys
import argparse
import os
import pysam

### def print_aln(aln):
###     print "the alignment range is 0-based, and [ , )"
###     print
###     print ">>> raw alignment data"
###     print aln
### 
###     print "=" * 80
###     print ">>> cigar string"
###     print aln.cigarstring
###     print ">>> mapping quality = ", aln.mapping_quality
### 
###     print "=" * 40 + " query " + "="*40
###     print ">>> query_name =", aln.query_name
###     print ">>> is reverse = ", aln.is_reverse
###     print ">>> flag = ", aln.flag
###     #print ">>> inferred read length from CIGAR string =", aln.infer_query_length(False)
###     print ">>> query length =", aln.query_length
###     print ">>> query length before partitioned into 3 parts =", aln.query_name.split('/', 1)[0].split('_')[-1]
###     print ">>> query alignment pos"
###     print "[{}, {}): {}".format(aln.query_alignment_start, aln.query_alignment_end, aln.query_alignment_length)
### 
###     print "=" * 40 + " reference " + "=" * 40
###     print ">>> reference id =", aln.reference_id
###     print ">>> reference alignment pos"
###     print "[{}, {}): {}".format(aln.reference_start, aln.reference_end, aln.reference_length)

def run(args):
    fin_bam = pysam.AlignmentFile(args.bam, "rb")
    fout = open(os.path.join(args.directory, "brief_alignment"), "w") if args.directory else sys.stdout
    # write total # of references
    fout.write("{}\n".format( fin_bam.nreferences ))
    # write reference name: <length> <reference name>
    lengths = fin_bam.lengths
    for i in xrange(fin_bam.nreferences):
        fout.write("{} {}\n".format( lengths[i], fin_bam.get_reference_name(i) ))
    # write simple alignment
    for aln in fin_bam.fetch(until_eof = True):
        qname = aln.query_name
        if args.modify_qname:
            qname = qname.split("/", 1)[0]
            qname_tokens = qname.split("_")
            qname = "{}/{}/0_{}".format(qname, qname_tokens[0][5:], int(qname_tokens[2])-1)
        fout.write("{} {} {} ".format(aln.reference_id, qname, aln.query_length))
        fout.write("{} {} {} {} {} {}\n".format(aln.reference_start, aln.reference_end, 
                aln.query_alignment_start, aln.query_alignment_end,
                aln.mapping_quality, 'R' if aln.is_reverse else 'F'))
    # write each alignment header
    if args.directory:
        fout.close()
    fin_bam.close()



def parse_args( argv = None ):
    parser = argparse.ArgumentParser(description = "Bam file extractor")
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    #parser.add_argument("-t", "--target", required=True, help="contig/scaffold file (fasta/q)")
    #parser.add_argument("-r", "--read", required=True, help="long reads file (fasta/q)")
    parser.add_argument("-d", "--directory", help="directory of the output file. If not set, it will be output to stdout")
    parser.add_argument("-m", "--modify-qname", action="store_true", help="modify query name")
    return parser.parse_args( argv )

def main( argv = None ):
    args = parse_args( argv )
    run( args )
    

if __name__ == "__main__":
    #main(["-b", "/Volumes/bioinfo/tmp/finding_inversions/afun_out.ctg.bam", "-d", ".", "-m"])
    main()
