
import os
import json
import sys
import traceback
import argparse
import common
from subprocess import call

## Please change the gtf filename here if you want to use an updated version.
reffile = "references/Mus_musculus.NCBIM37.63.chr.gtf"

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
    common.run('cuffcompare -r ' + reffile + ' -o ' +
               resultsubfolder  + "/cuffcompare " +
               cuffresultfolder + '/' + alias + '/transcripts.gtf')
