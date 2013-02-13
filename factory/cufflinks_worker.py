
import os
import json
import sys
import traceback
import argparse
import common
from subprocess import call

def cufflinks_worker(args):
    resultfolder = 'result/' + os.path.basename(args.filelist) + \
                   args.cufflinks_folder
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    f = args.filename
    resultsubfolder = resultfolder + '/' + \
                      os.path.basename(f.strip().replace('.bam',''))
    if not os.path.exists(resultsubfolder):
        os.makedirs(resultsubfolder)
    if (args.is_new_cufflinks):
        common.run('cufflinks -F 0.1 -p 1 -q -o ' + resultsubfolder +
                   ' ' + f.strip())
    else :
        common.run('old_cufflinks  -F 0.1 -p 1 -q  -o ' + resultsubfolder +
                   ' ' + f.strip())
