#!/usr/bin/env python

import os
import invdet
from invdet import bamExtractor
from invdet import peGenerator
from invdet.invdet_core import InvDector
from invdet.maxcut import MaxCut
import errno
import logging
import itertools
import argparse
import subprocess
import pysam
import multiprocessing
import signal

def nucmer_task(params, logger=None):
    log_fname = os.path.join( params[1], "nucmer.{}.log".format( params[2]) )
    prefix = os.path.join( params[1], "nc_aln.{}".format( params[2] ) )
    infile_name = os.path.join( params[1], "nucmer-input.{}.fasta".format( params[2] ) )
    with open(log_fname, "wb") as fout_nclog:
        child_process = subprocess.Popen( params[0] + ["--prefix", prefix, infile_name, infile_name], stderr = fout_nclog )
        ret_code = child_process.wait()
        if ret_code != 0:
            if logger:
                logger.error("Run Nucmer failed. See `{}` for more details".format( log_fname ))
                raise subprocess.CalledProcessError(ret_code, params[0] + ["--prefix", prefix, infile_name, infile_name])
            else:
                raise RuntimeError("[ERROR]: Run Nucmer failed. See `{}` for more details. ret_code = {}, command = {}".format( log_fname, ret_code, str(params[0] + ["--prefix", prefix, infile_name, infile_name])))

def makedir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def parse_args(argv = None):
    stage_list = ["begin", "blasr", "extract", "inv-repeats", "report"]
    parser = argparse.ArgumentParser(description = "Inversion Dector")
    parser.add_argument("-d", "--working-directory", default="working_dir", help="working directory (default: %(default)s)")
    parser.add_argument("-t", "--target-genome", help="Target genome file (fasta file)")
    parser.add_argument("-r", "--reads", help="reads file (fasta/fastq format, pacbio long reads preferred)")
    parser.add_argument("-s", "--start-from", default="begin", choices=stage_list, help="start the program from (default: %(default)s)")
    parser.add_argument("-o", "--only", action="store_true", help="Run the chosen stage only")
    parser.add_argument("-S", "--chosen-stages", action="append", choices=stage_list, help="Choose a stage or stages to run")
    parser.add_argument("--strategy", default="ignore-IR", choices=["naive", "ignore-IR", "extract-IR"], help="Choose a strategy (default: %(default)s) ('IR' stands for 'Inverted Repeats'; 'extract-IR' has not been implemented yet)")
    parser.add_argument("-V", "--version", action='version', version=('%(prog)s ' + invdet.__version__))
    parser.add_argument("--min-coverage", default=5, type=int, help="min coverage for filtering poor alignments (default: %(default)s)")
    parser.add_argument("--min-percent", default=0.05, type=float, help="min percentage of coverage for filtering poor alignments (default: %(default)s)")
    parser.add_argument("--min-overlap", default=50, type=int, help="min overlap for determining the overlaped regions (default: %(default)s)")
    parser.add_argument("--small-graph", default=15, type=int, help="max number of vertices for small graph (default: %(default)s)")
    parser.add_argument("--max-nodes", default=60, type=int, help="max number of vertices that can run on with 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--min-iter", default=100, type=int, help="min iterations for running 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--max-iter", default=10000, type=int, help="max iterations for running 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--min-ratio", default=0.878, type=float, help="min approx ratio for the 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--max-ratio", default=0.995, type=float, help="max approx ratio for the 0.878-approx algorithm (default: %(default)s)")
    parser.add_argument("--log", action="store_true", help="save log to file [invdet.log] instead of printing in the console")
    
    # blasr/nucmer options
    parser.add_argument("-j", "--nproc", default=1, type=int, help="BLASR/NUCmer option: number of threads (default: %(default)s)")

    # blasr options
    parser.add_argument("--minMatch", type=int, default=12, help="BLASR option: minimum seed length (default: %(default)s)")
    parser.add_argument("--minReadLength", type=int, default=50, help="BLASR option: Skip reads that have a full length less than minReadLength. Subreads may be shorter. (default: %(default)s)")
    parser.add_argument("--minAlnLength", type=int, default=0, help="BLASR option: report alignments only if their lengths are greater than minAlnLength (default: %(default)s)")
    parser.add_argument("--minPctSimilarity", type=float, default=0, help="BLASR option: report alignments only if their percentage similarity is greater than minPctSimilarity (range: [0.0, 100.0], default: %(default)s)")
    parser.add_argument("--minPctAccuracy", type=float, default=0, help="BLASR option: report alignments only if their percentage accuracy is greater than minPctAccuracy (range: [0.0, 100.0], default: %(default)s)")
    
    # nucmer options
    parser.add_argument("--nc-breaklen", type=int, default=50, help="Nucmer option: set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default: %(default)s)")
    parser.add_argument("--nc-mincluster", type=int, default=65, help="Nucmer option: set the minimum length of a cluster of matches (default: %(default)s)")
    parser.add_argument("--nc-diagdiff", type=int, default=5, help="Nucmer option: set the maximum diagonal difference between two adjacent anchors in a cluster (default: %(default)s)")
    parser.add_argument("--nc-diagfactor", type=float, default=0.12, help="Nucmer option: set the maximum diagonal difference between two adjacent anchors in a cluster as a differential fraction of the gap length (default: %(default)s)")
    parser.add_argument("--nc-maxgap", type=int, default=90, help="Nucmer option: set the maximum gap between two adjacent matches in a cluster (default: %(default)s)")
    parser.add_argument("--nc-minmatch", type=int, default=20, help="Nucmer option: Set the minimum length of a single match (default: %(default)s)")

    return parser.parse_args( argv )

def addCmdParameter(command, param, *parameters):
    command.append( param )
    for each in parameters:
        command.append( str(each) )

def run_begin(args, logger):
    logger.info("[begin] Generate PE files from reads")
    if not args.reads:
        logger.critical("-r/--reads is missing")
        exit(-1)
    peGenerator.main(["-i", args.reads, "-p", os.path.join( args.working_directory, "pe_reads" ), "-t", "fastx", "-f", "afun", "-s", "-P"])

def run_blasr(args, logger):
    logger.info("[blasr] Call Blasr to align the long reads to the target genome")
    if not args.target_genome:
        logger.critical("-t/--target-genome is missing")
        exit(-1)
    pe_prefix = os.path.join( args.working_directory, "pe_reads")
    blasr_command = ["blasr", pe_prefix + ".fastq", args.target_genome, "--allowAdjacentIndels", "--out", pe_prefix+".bam", "--bam", "--unaligned", pe_prefix + ".unaligned.txt", "--noPrintUnalignedSeqs", "--clipping", "hard", "--hitPolicy", "allbest"]
    if args.nproc > 1:
        addCmdParameter(blasr_command, "--nproc", args.nproc)
    if args.minMatch != 12:
        addCmdParameter(blasr_command, "--minMatch", args.minMatch)
    if args.minReadLength != 50:
        addCmdParameter(blasr_command, "--minReadLength", args.minReadLength)
    if args.minAlnLength > 0:
        addCmdParameter(blasr_command, "--minAlnLength", args.minAlnLength)
    if args.minPctSimilarity > 0:
        if args.minPctSimilarity > 100:
            logger.error("--minPctSimilarity should be in the range [0.0, 100.0]")
            exit(1)
        addCmdParameter(blasr_command, "--minPctSimilarity", args.minPctSimilarity)
    if args.minPctAccuracy > 0:
        if args.minPctAccuracy > 100:
            logger.error("--minPctAccuracy should be in the range [0.0, 100.0]")
            exit(1)
        addCmdParameter(blasr_command, "--minPctAccuracy", args.minPctAccuracy)
    logger.debug("blasr parameter: {}".format(str(blasr_command)))
    subprocess.check_call(blasr_command)

def run_inv_repeats(args, logger):
    if args.strategy == "naive":
        return
    logger.info("[inv-repeats] Extract inverted repeats using Nucmer")
    nucmer_working_dir = os.path.join( args.working_directory, "nucmer" )
    makedir( nucmer_working_dir )
    if not args.target_genome:
        logger.critical("-t/--target-genome is missing")
        exit(-1)
    nucmer_command = ["nucmer", "--maxmatch", "--reverse"]
    if args.nc_breaklen != 200:
        addCmdParameter(nucmer_command, "-b", args.nc_breaklen)
    if args.nc_mincluster != 65:
        addCmdParameter(nucmer_command, "-c", args.nc_mincluster)
    if args.nc_diagdiff != 5:
        addCmdParameter(nucmer_command, "-D", args.nc_diagdiff)
    if args.nc_maxgap != 90:
        addCmdParameter(nucmer_command, "-g", args.nc_maxgap)
    if args.nc_minmatch != 20:
        addCmdParameter(nucmer_command, "-l", args.nc_minmatch)
    #addCmdParameter(nucmer_command, args.target_genome, args.target_genome)
    logger.debug("numcer parameters: {}".format(str(nucmer_command)))

    n_refs = 0
    with pysam.FastxFile(args.target_genome) as fh:
        nucmer_input_prefix = os.path.join( nucmer_working_dir, "nucmer-input" )
        for entry in fh:
            with open(nucmer_input_prefix + ".{}.fasta".format(n_refs), "w") as fout:
                fout.write(">{}\n".format( n_refs ))
                fout.write(entry.sequence)
                fout.write("\n")
            n_refs += 1
    with open(os.path.join( nucmer_working_dir, "n_refs.txt"), "w") as fout:
        fout.write("{}\n".format( n_refs ))
    nthreads = min(n_refs, args.nproc)
    if nthreads <= 1:
        logger.info("run nucmer with single thread on {} targets".format( n_refs ))
        for i in xrange( n_refs ):
            nucmer_task( [nucmer_command, nucmer_working_dir, i], logger )
    else:
        task_params = [(nucmer_command, nucmer_working_dir, i) for i in xrange(n_refs)]
        original_sigint_handler = signal.signal( signal.SIGINT, signal.SIG_IGN )
        pool = multiprocessing.Pool( nthreads )
        signal.signal( signal.SIGINT, original_sigint_handler )
        
        logger.info("run nucmer with {} threads on {} targets, with blocksize={}".format(nthreads, n_refs, min(128, n_refs/nthreads)))
        res = pool.map_async(nucmer_task, task_params, min(128, n_refs/nthreads))
        while True:
            try:
                res.get(0x7fffffff)
            except KeyboardInterrupt:
                pool.terminate()
                pool.join()
                raise
            except multiprocessing.TimeoutError:
                pass
            except RuntimeError:
                pool.terminate()
                pool.join()
                raise
            else:
                pool.close()
                break
        pool.join()
    
def run_extract(args, logger): 
    logger.info("[extract] Extract alignments from BAM file")
    bamExtractor.main(["-b", os.path.join(args.working_directory, "pe_reads.bam"), "-d", args.working_directory])

def run_report(args, logger):
    logger.info("[report] Generate report")
    logger.info("generate graph")
    inv_dector = InvDector()
    brief_alignment = os.path.join( args.working_directory, "brief_alignment")
    graph_file = os.path.join( args.working_directory, "graph_file")
    inv_dector.read( brief_alignment )
    if args.strategy == "ignore-IR":
        inv_dector.gen_graphs_ignore_inverted_repeats(graph_file, os.path.join(args.working_directory, "nucmer"), args.min_coverage, args.min_percent, args.min_overlap)
    elif args.strategy == "naive":
        inv_dector.gen_graphs(graph_file, args.min_coverage, args.min_percent, args.min_overlap)
    else: # "extract-IR"
        logger.error("The strategy 'extract-IR' has not been implemented yet! Please choose another strategy")
        exit(1)

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
    

def main(argv = None):
    args = parse_args(argv)
    makedir( args.working_directory )

    if args.log:
        logging.basicConfig(filename=os.path.join(args.working_directory, "invdet.log"), format="[%(asctime)s] [%(levelname)s] %(message)s", level=logging.DEBUG)
    else:
        logging.basicConfig(format="[%(asctime)s] [%(levelname)s] %(message)s", level=logging.DEBUG)
    logger = logging.getLogger()
    logger.info("Start")

    func_list = [run_begin, run_blasr, run_extract, run_inv_repeats, run_report]
    func_indices = {"begin": 0, 
                    "blasr": 1, 
                    "extract": 2, 
                    "inv-repeats": 3, 
                    "report": 4}
    if args.only:
        stage_ids = set()
        for each in args.chosen_stages:
            stage_ids.add( func_indices[each] )
        stage_ids = list( stage_ids )
        stage_ids.sort()
        for each in stage_ids:
            func_list[ each ](args, logger)
    else:
        for each_func in func_list[ func_indices[args.start_from]: ]:
            each_func(args, logger)
    logger.info("[Done!]")

if __name__ == "__main__":
    main()

