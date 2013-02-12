import os
import json
import sys
import traceback
import argparse
from subprocess import call

def jeweler_worker(args, refidtable, reffiletable):
    files=open(args.filelist).readlines()
    result_folder = ""
    appraiser_result_folder = ""
    if (args.is_new_cufflinks):
        result_folder = 'result/'+os.path.basename(args.filelist)+'/new_jeweler/'
    else:
        result_folder = 'result/'+os.path.basename(args.filelist)+'/jeweler/'
    if (args.is_new_cufflinks):
        appraiser_result_folder = 'result/'+os.path.basename(args.filelist)+'/new_appraiser/'
    else:
        appraiser_result_folder = 'result/'+os.path.basename(args.filelist)+'/appraiser/'
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    for f in files:
        resultsubfolder=result_folder+'/'+os.path.basename(f.strip().replace('.bam',''))
        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)
        alias=os.path.basename(f.strip().replace('.bam',''))
        unmapped_left_fastq_file = f.strip().replace('_merged.bam', '.left.unmapped.fastq')
        unmapped_right_fastq_file = f.strip().replace('_merged.bam', '.right.unmapped.fastq')
        maternal_strain_id=refidtable[alias[0]]
        paternal_strain_id=refidtable[alias[1]]
        maternal_strain_ref_seq=reffiletable[alias[0]]
        paternal_strain_ref_seq=reffiletable[alias[1]]
        bam_file=f.strip()

        result_file=alias+'.info'
        cuffresult_folder='result/'+os.path.basename(args.filelist)+args.cufflinks_folder
        gtf_input_file=	cuffresult_folder+'/'+os.path.basename(f.strip().replace('.bam',''))+'/transcripts.gtf'
        jeweler_args = (' -maternal_strain_ref_file '+maternal_strain_ref_seq +' ' +
                        '-paternal_strain_ref_file '+paternal_strain_ref_seq +' ' +
                        '-maternal_strain_id '+maternal_strain_id +' ' +
                        '-paternal_strain_id '+paternal_strain_id +' ' +
                        '-bam_file '+bam_file +' ' +
                        '-alias '+alias +' ' +
                        '-result_file '+result_file +' ' +
                        '-result_folder '+result_folder +' '+
                        '-gtf_input_file '+gtf_input_file+' ' +
                        "-left_unmapped_file " + unmapped_left_fastq_file + " "+
                        "-right_unmapped_file " + unmapped_right_fastq_file  )
        jeweler_command = './jeweler -i '+result_folder+'/'+alias+'/'+result_file
        mamf_command = ""
        bracelet_command = ""
        earrings_command = ""
        ma_command = ""
        prepare_command = ""
        # if (os.path.exists(appraiser_result_folder + alias +"/sm_map")):
        mamf_command = " -mamf " + appraiser_result_folder + alias +"/sm_map";

        if (args.is_earrings):
            ## multiple alignments maps
            earrings_command = " -earrings"
        if (args.is_prepare):
            ## multiple alignments maps
            prepare_command = " -prepare"
        if (args.is_mismatch_analyzer):
            if (args.is_new_cufflinks):
                ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/new_mismatch_analyzer/' + alias + '/'
            else:
                ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/' + alias + '/'
            if not os.path.exists(ma_result_folder):
                os.makedirs(ma_result_folder)
            ma_command = " -mismatch_analyzer " + ma_result_folder+ "result"
        if (args.is_bracelet):
            ma_result_folder = ""
            if (args.is_new_cufflinks):
                ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/new_bracelet/' + alias + '/'
            else:
                ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/bracelet/' + alias + '/'
            if not os.path.exists(ma_result_folder):
                    os.makedirs(ma_result_folder)
            bracelet_command = " -bracelet " + ma_result_folder+ " "
        jeweler_command = jeweler_command + mamf_command + earrings_command + bracelet_command+ ma_command+ jeweler_args + prepare_command
        print(jeweler_command)