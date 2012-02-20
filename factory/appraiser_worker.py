import os
import json
import sys
import traceback
import argparse
from subprocess import call

def appraiser_worker(args):
    files=open(args.filelist).readlines()
    resultfolder='result/'+os.path.basename(args.filelist)+'/appraiser/'
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    for f in files:
        resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)
        print('./appraiser -mamf '+resultsubfolder+'/sm -l '+resultsubfolder+'/log -qf '+ resultsubfolder +'/qf -b '+f.strip())

