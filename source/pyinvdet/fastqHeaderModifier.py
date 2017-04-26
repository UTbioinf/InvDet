#!/usr/bin/env python

import argparse

def parse_args( argv = None):
    parser = argparse.ArgumentParser(description="change header of fastq files")
    parser.add_argument("infile", help="input fastq file name")
    parser.add_argument("outfile", help="output fastq file name")
    return parser.parse_args( argv )

def main( argv = None ):
    args = parse_args()
    fin = open(args.infile, "r")
    fout = open(args.outfile, "w")
    fout_log = open(args.outfile + ".log", "w")
    while True:
        line = fin.readline()
        if not line:    break
        # generate new header
        line = line.strip()
        line_tokens = line.split('_')
        qid = line_tokens[0][line_tokens[0].find("afun")+len("afun"):]
        new_header = "{}/{}/0_{}\n".format( line, qid, int(line_tokens[2]) - 1 )
        fout.write(new_header)
        fout_log.write("{}\n{}\n".format(line, new_header))

        # write sequence
        line = fin.readline()
        fout.write( line )
        # write +[tag]
        line = fin.readline()
        if line[1] == '\n' or line[1] == '\r':
            fout.write( line )
        else:
            fout.write("+{}\n".format(new_header))
        # write quality
        line = fin.readline()
        fout.write( line )
    fin.close()
    fout.close()
    fout_log.close()


if __name__ == "__main__":
    main()
