
import os
import json
import sys
import traceback
import argparse
from subprocess import call

def cuffcompare_worker(args):
    files=open(args.filelist).readlines()
    ## run cuffcompare command
    resultfolder='result/'+os.path.basename(args.filelist)+'/cuffcompare'
    cuffresultfolder='result/'+os.path.basename(args.filelist)+'/cufflinks/'
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    for f in files:
        resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))

        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)

        alias=os.path.basename(f.strip().replace('.bam',''))
        print('cuffcompare -r data/database/Mus_musculus.NCBIM37.63.chr.gtf -o '+resultsubfolder + '/cuffcompare '+cuffresultfolder+'/'+alias+'/transcripts.gtf')


                