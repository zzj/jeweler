import os
import argparse
import sys
import traceback


def initialize_parser():
    parser=argparse.ArgumentParser(description=
                                   'The manager for RNA-seq analysis.')

    parser.add_argument('--filelist',
                        help='A list of bam files, for multipe bamfiles. ',
                        dest='filelist')
    
    parser.add_argument('--reftable',
                        help='A reference table for reference genomes',
                        dest='reftable')

    parser.add_argument('--new_cufflinks',
                        action='store_true',
                        dest='is_new_cufflinks')

    return parser


def run_command(file_list, reftable, parameters, extra):
    cmd = "python3.2 factory/miner.py --filelist " + file_list
    cmd += " --reftable " + reftable
    cmd += " " + parameters + " " + extra
    cmd += ">> command_list"
    print(cmd)
    os.system(cmd)


def pause_command():
    os.system("echo pause >> command_list")
    

def main():
    try:
        parser = initialize_parser()
        args = parser.parse_args();
        if not args.filelist or not args.reftable:
            parser.error("Require filelist and reftable")
        extra = ""
        if args.is_new_cufflinks:
            extra = "--new_cufflinks"
        os.system("rm command_list")
        run_command(args.filelist, args.reftable, "--cufflinks", extra)
        run_command(args.filelist, args.reftable, "--appraiser", extra)
        pause_command()
        run_command(args.filelist, args.reftable, "--cuffcompare", extra)
        run_command(args.filelist, args.reftable,
                    "--jeweler --earrings --bracelet", extra)
        pause_command()
        run_command(args.filelist, args.reftable,
                    "--jeweler --mismatch_analyzer", extra)
        pause_command()
        run_command(args.filelist, args.reftable, "--shared_graph", extra)
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("*** print_exc:")
        traceback.print_exc()


if __name__ == "__main__":
    main()
