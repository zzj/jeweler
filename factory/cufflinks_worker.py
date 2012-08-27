
import os
import json
import sys
import traceback
import argparse
from subprocess import call

def cufflinks_worker(args):
                ## run cufflinks command
    files=open(args.filelist).readlines()
    resultfolder='result/'+os.path.basename(args.filelist)+args.cufflinks_folder
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    for f in files:
        resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)
        if (args.is_new_cufflinks):
            print('cufflinks -F 0.1 -p 1 -q -o '+resultsubfolder+' '+f.strip())
        else :
            print('old_cufflinks  -F 0.1 -p 1 -q  -o '+resultsubfolder+' '+f.strip())
