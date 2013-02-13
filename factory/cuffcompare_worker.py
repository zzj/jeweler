
import os
import json
import sys
import traceback
import argparse
import common
from subprocess import call

def cuffcompare_worker(args):
    ## run cuffcompare command
    resultfolder='result/' + os.path.basename(args.filelist) + args.cuffcompare_folder
    cuffresultfolder='result/'+ os.path.basename(args.filelist) + args.cufflinks_folder
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    f = args.filename
    resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))

    if not os.path.exists(resultsubfolder):
        os.makedirs(resultsubfolder)

    alias=os.path.basename(f.strip().replace('.bam',''))
    common.run('cuffcompare -r data/database/Mus_musculus.NCBIM37.63.chr.gtf -o ' +
               resultsubfolder  + "/cuffcompare " +
               cuffresultfolder + '/' + alias + '/transcripts.gtf')
