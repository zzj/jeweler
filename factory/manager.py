import os
import argparse
import sys
import traceback

def config(queue, memory, process):
    print("config " + queue + " " + memory + " " + process)


def initialize_parser():
    parser=argparse.ArgumentParser(description=
                                   'The manager for RNA-seq analysis.')

    parser.add_argument('--filelist',
                        help='A bam file path. ',
                        dest='filelist')

    parser.add_argument('--id',
                        help='The line index in the filelist. ',
                        dest='id')

    parser.add_argument('--reftable',
                        help='A reference table for reference genomes',
                        dest='reftable')

    parser.add_argument('--new_cufflinks',
                        action='store_true',
                        dest='is_new_cufflinks')

    return parser


def run_command(filelist, filename, reftable, parameters, extra):
    cmd = "python3.2 factory/miner.py --filelist " + filelist + \
          " --reftable " + reftable + ' --filename ' + filename + \
          " " + parameters + " " + extra
    os.system(cmd)

def run_manager_command(filelist, i, reftable, extra):
    cmd = "python3.2 factory/manager.py --filelist " + filelist + \
          " --reftable " + reftable + ' --id ' + str(i) + ' ' + extra
    print(cmd)

def main():
    try:
        parser = initialize_parser()
        args = parser.parse_args();
        if not args.filelist or not args.reftable:
            parser.error("Require --filelist and --reftable and --id")
        extra = ""
        files = open(args.filelist, 'r').readlines()
        if args.is_new_cufflinks:
            extra = "--new_cufflinks"
        if args.id is not None:
            jobid = int(args.id)
            filename = files[jobid].strip()
            # run_command(args.filelist, filename, args.reftable,
            #             "--cufflinks", extra)
            # run_command(args.filelist, filename, args.reftable,
            #             "--appraiser", extra)
            # run_command(args.filelist, filename, args.reftable,
            #             "--cuffcompare", extra)
            # run_command(args.filelist, filename, args.reftable,
            #             "--jeweler --earrings --bracelet", extra)
            # run_command(args.filelist, filename, args.reftable,
            #             "--jeweler --mismatch_analyzer", extra)
            # run_command(args.filelist, filename, args.reftable,
            #             "--shared_graph", extra)
            run_command(args.filelist, filename, args.reftable,
                        "--classify_gene", extra)
        else:
            config("week", "12", "1")
            for i in range(len(files)):
                run_manager_command(args.filelist, i, args.reftable, extra)
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("*** print_exc:")
        traceback.print_exc()


if __name__ == "__main__":
    main()
