
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
        ## Sometimes, when the jobs scheduled at the same time, they
        ## will try to create folder at the same time, therefore, we
        ## need this "mkdir -p" function to avoid the exception.
        common.mkdir_p(resultfolder)
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
