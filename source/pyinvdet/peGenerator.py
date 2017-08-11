#!/usr/bin/env python

import pysam
import math
import argparse

class FileWriter(object):
    def __init__(self, write_type = "two-file", prefix = "output", use_pacbio_head = False):
        self._write_type = write_type
        self._prefix = prefix
        self._fout1 = None
        self._fout2 = None
        self._index = 0
        self._write = self._get_write_func( write_type )
        self._use_pacbio_head = use_pacbio_head

    def set_write_type(self, write_type):
        self._write_type = write_type
        self._write = self._get_write_func( write_type )

    def set_prefix(self, prefix):
        self._prefix = preifx


    def open(self):
        if self._write_type == "two-file" or self._write_type == "MP":
            self._fout1 = open(self._prefix + "_1.fastq", "w")
            self._fout2 = open(self._prefix + "_2.fastq", "w")
        else:
            self._fout1 = open(self._prefix + ".fastq", "w")
        self._index = 0

    def close(self):
        if self._fout1:
            self._fout1.close()
            self._fout1 = None
        if self._fout2:
            self._fout2.close()
            self._fout2 = None

    def write(self, head, head_remaining, seq, qual = None):
        self._write(head, ' ' + head_remaining if head_remaining else "", seq, qual)

    def _get_write_func(self, write_type):
        if write_type == "two-file":
            return self._write_two_file
        elif write_type == "single-file":
            return self._write_single_file
        elif write_type == "afun":
            return self._write_afun
        elif write_type == "PE":
            return self._write_PE
        elif write_type == "MP":
            return self._write_MP
        else:
            raise ValueError("Unknown file type '{}'to write".format( write_type ))

    def _nuc_rc(self, ch):
        if ch == 'A': return 'T'
        if ch == 'a': return 't'
        if ch == 'T': return 'A'
        if ch == 't': return 'a'
        if ch == 'G': return 'C'
        if ch == 'g': return 'c'
        if ch == 'C': return 'G'
        if ch == 'c': return 'g'
        if ch == 'N': return 'N'
        if ch == 'n': return 'n'
        raise ValueError("unknown base [{}]".format(ch))
    
    def _gen_rc(self, seq):
        return "".join([ self._nuc_rc( each ) for each in seq[::-1] ])

    def _write_two_file(self, head, head_remaining, seq, qual = None):
        length = len(seq)/3
        if length == 0: return
        self._write_fastq("@{}/1{}\n".format(head, head_remaining),
                seq[:length], qual[:length] if qual else ("!" * length), self._fout1)

        self._write_fastq("@{}/2{}\n".format(head, head_remaining),
                seq[-length:], qual[-length:] if qual else ("!" * length), self._fout2)

    def _write_MP(self, head, head_remaining, seq, qual = None):
        length = len(seq)/3
        if length == 0: return
        self._write_fastq("@{}/1{}\n".format(head, head_remaining),
                seq[:length], qual[:length] if qual else ("!" * length), self._fout1)

        self._write_fastq("@{}/2{}\n".format(head, head_remaining),
                self._gen_rc( seq[-length:]) , qual[:-length-1:-1] if qual else ("!" * length), self._fout2)

    def _write_single_file(self, head, head_remaining, seq, qual = None):
        length = len(seq) / 3
        if length == 0: return
        self._write_fastq("@{}/1{}\n".format(head, head_remaining),
                seq[:length], qual[:length] if qual else ("!" * length), self._fout1)

        self._write_fastq("@{}/2{}\n".format(head, head_remaining),
                seq[-length:], qual[-length:] if qual else ("!" * length), self._fout1)

    def _write_afun(self, head, head_remaining, seq, qual = None):
        total_len = len(seq)
        length = total_len / 3
        if length == 0: return
        pacbio_header = "/{}/0_{}".format( self._index, length ) if self._use_pacbio_head else ""
        self._write_fastq("@afun{}_5_{}{}\n".format(self._index, total_len, pacbio_header), 
                seq[:length], qual[:length] if qual else ("!" * length), self._fout1)

        pacbio_header = "/{}/{}_{}".format( self._index, total_len - length, total_len ) if self._use_pacbio_head else ""
        self._write_fastq("@afun{}_3_{}{}\n".format(self._index, total_len, pacbio_header),
                seq[-length:], qual[-length:] if qual else ("!" * length), self._fout1)
        self._index += 1

    def _write_PE(self, head, head_remaining, seq, qual = None):
        total_len = len(seq)
        length = total_len / 3
        if length == 0: return
        pacbio_header = "/{}/0_{}".format( self._index, length ) if self._use_pacbio_head else ""
        self._write_fastq("@PE_{}_5_{}{}\n".format(self._index, total_len, pacbio_header), 
                seq[:length], qual[:length] if qual else ("!" * length), self._fout1)

        pacbio_header = "/{}/{}_{}".format( self._index, total_len - length, total_len ) if self._use_pacbio_head else ""
        self._write_fastq("@PE_{}_3_{}{}\n".format(self._index, total_len, pacbio_header),
                seq[-length:], qual[-length:] if qual else ("!" * length), self._fout1)
        self._index += 1



    def _write_fastq(self, head_with_newline, seq, qual, fout):
        fout.write(head_with_newline)
        fout.write(seq)
        fout.write("\n+\n")
        fout.write(qual)
        fout.write("\n")

class Statistics(object):
    def __init__(self):
        self._cnt = 0
        self._len = 0.0
        self._square_len = 0.0
        self._gap_len = 0.0
        self._square_gap_len = 0.0

    def add_entry(self, entry):
        self._cnt += 1
        length = len(entry.sequence)
        self._len += length
        self._square_len += length * length
        length -= length / 3 * 2
        self._gap_len += length
        self._square_gap_len += length * length

    def write(self, fout):
        fout.write("[raw reads]\n")
        fout.write("count:              {}\n".format(self._cnt))
        fout.write("total length:       {}\n".format(int(self._len)))
        fout.write("average length:     {}\n".format(self._len / self._cnt))
        variance = self._square_len / self._cnt - self._len * self._len / (self._cnt * self._cnt)
        fout.write("variance:           {}\n".format(variance))
        fout.write("standard deviation: {}\n".format( math.sqrt(variance) ))

        fout.write("[gap]\n")
        fout.write("total length:       {}\n".format(int(self._gap_len)))
        fout.write("average length:     {}\n".format(self._gap_len / self._cnt))
        variance = self._square_gap_len / self._cnt - self._gap_len * self._gap_len / (self._cnt * self._cnt)
        fout.write("variance:           {}\n".format(variance))
        fout.write("standard deviation: {}\n".format( math.sqrt(variance) ))

def parse_args( argv = None):
    parser = argparse.ArgumentParser(description = "Generate PE reads from fasta/fastq file")
    parser.add_argument("-i", "--input-file", required=True, help="Input file name")
    parser.add_argument("-p", "--prefix", default="output", help="Prefix for the output files")
    parser.add_argument("-t", "--file-type", default="fastx", choices=["fastx"], help="Input file type (fastx includes both fasta and fastq, and probably a mixture of them)")
    parser.add_argument("-f", "--output-format", default="two-file", choices=["two-file", "single-file", "afun", "PE", "MP"], help="Output format)")
    parser.add_argument("-s", "--statistics", action="store_true", help="Generate statistics")
    parser.add_argument("-P", "--pacbio-head", action="store_true", help="Use pacbio head, which is only useful for `afun` and `PE` format")
    return parser.parse_args(argv)

def main(argv = None):
    args = parse_args( argv )

    fout = FileWriter(args.output_format, args.prefix, args.pacbio_head)
    fout.open()

    if args.statistics:
        stats = Statistics()
    with pysam.FastxFile(args.input_file) as fh:
        for entry in fh:
            fout.write(entry.name, entry.comment, entry.sequence, entry.quality)
            if args.statistics:
                stats.add_entry( entry )
    fout.close()
    if args.statistics:
        with open(args.prefix + ".stats.txt", "w") as fout:
            stats.write( fout )

if __name__ == "__main__":
    main()
