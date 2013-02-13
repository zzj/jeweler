import os
import json
import sys
import traceback
import argparse
import common
from subprocess import call

def appraiser_worker(args):

    if (args.is_new_cufflinks):
        resultfolder = 'result/' + os.path.basename(args.filelist) + '/new_appraiser/'
    else:
        resultfolder = 'result/' + os.path.basename(args.filelist) + '/appraiser/'
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)

    f = args.filename.strip()
    resultsubfolder=resultfolder+'/'+os.path.basename(f.replace('.bam',''))
    if not os.path.exists(resultsubfolder):
        os.makedirs(resultsubfolder)
    common.run('./appraiser -mamf ' + resultsubfolder + '/sm ' +
          '-l ' + resultsubfolder + '/log ' +
          '-qf '+ resultsubfolder +'/qf -b '+ f.strip())
