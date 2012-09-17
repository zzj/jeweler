import os
import json
import sys
import traceback
import argparse
from subprocess import call

from shared_graph_worker import shared_graph_worker
from shared_graph_worker import classify_gene_worker
from miner_parser import initialize_parser
from cufflinks_worker import cufflinks_worker
from cuffcompare_worker import cuffcompare_worker
from jeweler_worker import jeweler_worker
from appraiser_worker import appraiser_worker
from transcriptome_alignment_worker import transcriptome_alignment_worker

def config(queue, memory, process):
    print("config " + queue + " " + memory + " " + process)

def main():
    try:
        parser = initialize_parser()

        args=parser.parse_args();
        refidtable=dict()
        reffiletable=dict()
        if (args.reftable is not None and os.path.exists(args.reftable)):
            with open(args.reftable) as f:
                refs = f.readlines()
            for ref in refs:
                r=ref.strip().split('\t')
                refidtable[r[0]]=r[1]
                reffiletable[r[0]]=r[2]

        else :
            print('No reftable supplied')
            return

        if (args.filelist != None and os.path.exists(args.filelist)):
            # if (args.reftable == None):
            # 	print('Error: no reftable supplied (--reftable)')
            # 	return
            if (args.is_new_cufflinks):
                args.cufflinks_folder="/new_cufflinks/"
                args.cuffcompare_folder="/new_cuffcompare/"
                args.jeweler_folder="/new_jeweler/"
                args.bracelet_folder="/new_bracelet/"
                args.shared_graph_folder="/new_shared_graph/"
            else:
                args.cufflinks_folder="/cufflinks/"
                args.jeweler_folder="/jeweler/"
                args.cuffcompare_folder="/cuffcompare/"
                args.bracelet_folder="/bracelet/"
                args.shared_graph_folder="/shared_graph/"

            if (args.is_cufflinks):
                config("day", "8", "1")
                cufflinks_worker(args)
            elif (args.is_cuffcompare):
                config("day", "8", "1")
                cuffcompare_worker(args)
            elif (args.is_jeweler):
                config("day", "24", "1")
                jeweler_worker(args, refidtable, reffiletable)
            elif (args.is_transcriptome_alignment):
                transcriptome_alignment_worker(args, refidtable, reffiletable)
            elif (args.is_appraiser):
                config("day", "4", "1")
                appraiser_worker(args)
            elif (args.plot_shared_graph):
                config("day", "4", "1")
                shared_graph_worker(args)
            elif (args.classify_gene):
                classify_gene_worker(args)

    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("*** print_exc:")
        traceback.print_exc()

if __name__ == "__main__":
	main()
