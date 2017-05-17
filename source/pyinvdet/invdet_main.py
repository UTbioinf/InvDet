#!/usr/bin/env python

import os
from invdet import bamExtractor
from invdet import peGenerator
from invdet.invdet_core import InvDector
from invdet.maxcut import MaxCut
import errno
import logging
import itertools
import argparse
import subprocess

def makedir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def parse_args(argv = None):
    parser = argparse.ArgumentParser(description = "Inversion Dector")
    parser.add_argument("-d", "--working-directory", default="working_dir", help="working directory (default: %(default)s)")
    parser.add_argument("-t", "--target-genome", help="Target genome file")
    parser.add_argument("-r", "--reads", help="reads file (fasta/fastq format, pacbio long reads preferred)")
    parser.add_argument("-s", "--start-from", default="begin", choices=["begin", "blasr", "extract", "report"], help="start the program from (default: %(default)s)")
    parser.add_argument("-j", "--nproc", default=1, type=int, help="number of threads for BLASR (default: %(default)s)")
    parser.add_argument("--min-coverage", default=5, type=int, help="min coverage for filtering poor alignments (default: %(default)s)")
    parser.add_argument("--min-percent", default=0.05, type=float, help="min percentage of coverage for filtering poor alignments (default: %(default)s)")
    parser.add_argument("--small-graph", default=15, type=int, help="max number of vertices for small graph (default: %(default)s)")
    parser.add_argument("--max-nodes", default=60, type=int, help="max number of vertices that can run on with 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--min-iter", default=100, type=int, help="min iterations for running 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--max-iter", default=10000, type=int, help="max iterations for running 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--min-ratio", default=0.878, type=float, help="min approx ratio for the 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--max-ratio", default=0.995, type=float, help="max approx ratio for the 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--log", action="store_true", help="save log to file [invdet.log] instead of printing in the console")

    return parser.parse_args( argv )

def main(argv = None):
    args = parse_args(argv)
    makedir( args.working_directory )

    if args.log:
        logging.basicConfig(filename=os.path.join(args.working_directory, "invdet.log"), format="[%(asctime)s] [%(levelname)s] %(message)s", level=logging.DEBUG)
    else:
        logging.basicConfig(format="[%(asctime)s] [%(levelname)s] %(message)s", level=logging.DEBUG)
    logger = logging.getLogger()
    logger.info("Start")

    start_from = args.start_from
    if start_from == "begin":
        if not args.reads:
            logger.critical("-r/--reads is missing")
            exit(-1)
        logger.info("Generate PE files from reads")
        pe_prefix = os.path.join( args.working_directory, "pe_reads" )
        peGenerator.main(["-i", args.reads, "-p", pe_prefix, "-t", "fastx", "-f", "afun", "-s", "-P"])
        start_from = "blasr"

    if start_from == "blasr":
        if not args.target_genome:
            logger.critical("-t/--target-genome is missing")
            exit(-1)
        logger.info("Run BLASR")
        pe_prefix = os.path.join( args.working_directory, "pe_reads")
        subprocess.call(["blasr", pe_prefix + ".fastq", args.target_genome, "--allowAdjacentIndels", "--out", pe_prefix+".bam", "--bam", "--unaligned", pe_prefix + ".unaligned.txt", "--noPrintUnalignedSeqs", "--clipping", "hard", "--nproc", str(args.nproc)])
        start_from = "extract"

    if start_from == "extract":
        logger.info("extract bam")
        bamExtractor.main(["-b", os.path.join(args.working_directory, "pe_reads.bam"), "-d", args.working_directory])
        start_from = "report"

    if start_from == "report":
        logger.info("generate graphs")
        inv_dector = InvDector()
        brief_alignment = os.path.join( args.working_directory, "brief_alignment")
        graph_file = os.path.join( args.working_directory, "graph_file")
        inv_dector.read( brief_alignment )
        inv_dector.gen_graphs(graph_file, args.min_coverage, args.min_percent)

        logger.info("run max-cut")
        graph_cut = os.path.join(args.working_directory, "graph_cut")
        fout = open(graph_cut, "w")
        with open(graph_file, "r") as fin:
            while True:
                line = fin.readline()
                if not line: break
                line = line.split()
                n, r_id = int(line[0]), line[1]
                graph = MaxCut(args.small_graph, args.max_nodes)
                for i in xrange(n):
                    line = fin.readline().split()
                    graph.add_edge(int(line[0]), int(line[1]), int(line[2]))
                graph.solve(args.min_iter, args.max_iter, args.min_ratio, args.max_ratio)
                fout.write("{} {}\n".format(graph.number_of_nodes(), r_id))
                    
                for node_name, sol in itertools.izip(graph.node_name_iter(), graph.get_solution()):
                    fout.write("{} {}\n".format(node_name, 1 if sol else 0))
        fout.close()

        logger.info("Deduce inversions")
        inv_dector.report_inversions(graph_file, graph_cut,
                os.path.join(args.working_directory, "inversion.report"));
    logger.info("Done!")

if __name__ == "__main__":
    main()

