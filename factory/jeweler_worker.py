import os
import json
import sys
import traceback
import argparse
from subprocess import call

def jeweler_worker(args, refidtable, reffiletable):
    files=open(args.filelist).readlines()
    resultfolder = 'result/'+os.path.basename(args.filelist)+'/jeweler/'
    appraiser_resultfolder = 'result/'+os.path.basename(args.filelist)+'/appraiser/'
    if not os.path.exists(resultfolder):
        os.makedirs(resultfolder)
    for f in files:
        resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)
        result_folder=resultfolder
        alias=os.path.basename(f.strip().replace('.bam',''))
        maternal_strain_id=refidtable[alias[0]]
        paternal_strain_id=refidtable[alias[1]]
        maternal_strain_ref_seq=reffiletable[alias[0]]
        paternal_strain_ref_seq=reffiletable[alias[1]]
        bam_file=f.strip()
        result_file=alias+'.info'
        cuffresultfolder='result/'+os.path.basename(args.filelist)+args.cufflinks_folder
        gtf_input_file=	cuffresultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))+'/transcripts.gtf'
        jeweler_args = (' -maternal_strain_ref_file '+maternal_strain_ref_seq +' ' +
                        '-paternal_strain_ref_file '+paternal_strain_ref_seq +' ' +
                        '-maternal_strain_id '+maternal_strain_id +' ' +
                        '-paternal_strain_id '+paternal_strain_id +' ' +
                        '-bam_file '+bam_file +' ' +
                        '-alias '+alias +' ' +
                        '-result_file '+result_file +' ' +
                        '-result_folder '+result_folder +' '+
                        '-gtf_input_file '+gtf_input_file+' ')
        jeweler_command = './jeweler -i '+result_folder+'/'+alias+'/'+result_file
        mamf_command = ""
        bracelet_command = ""
        earrings_command = ""
        ma_command = ""
        if (os.path.exists(appraiser_resultfolder + alias +"/sm_table")):
            mamf_command = " -mamf " + appraiser_resultfolder + alias +"/sm_table";

        if (args.is_earrings):
            ## multiple alignments maps
            earrings_command = " -earrings" 
        if (args.is_mismatch_analyzer):
            ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/'
            if not os.path.exists(ma_resultfolder):
                os.makedirs(ma_resultfolder)
            ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/' + alias + '/'
            if not os.path.exists(ma_resultfolder):
                os.makedirs(ma_resultfolder)
            ma_command = " -mismatch_analyzer " + ma_resultfolder+ "result"
        if (args.is_bracelet):
            ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/bracelet/'
            if not os.path.exists(ma_resultfolder):
                    os.makedirs(ma_resultfolder)
            ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/bracelet/' + alias + '/'
            if not os.path.exists(ma_resultfolder):
                    os.makedirs(ma_resultfolder)
            bracelet_command = " -bracelet " + ma_resultfolder+ "result "
        jeweler_command = jeweler_command + mamf_command + earrings_command + bracelet_command+ ma_command+ jeweler_args
        print(jeweler_command)