import os
import json
import sys
import traceback
import argparse
from subprocess import call


from  miner_parser import initialize_parser
from cufflinks_worker import cufflinks_worker
from cuffcompare_worker import cuffcompare_worker
from jeweler_worker import jeweler_worker
from appraiser_worker import appraiser_worker

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
            

            if (args.is_cufflinks):
                cufflinks_worker(args)
            if (args.is_cuffcompare):
                cuffcompare_worker(args)
            elif (args.is_jeweler):
                jeweler_worker(args, refidtable, reffiletable)
            elif (args.is_appraiser) :
                appraiser_worker(args)
                
    except: 
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("*** print_exc:")
        traceback.print_exc()

if __name__ == "__main__":
	main()
		
