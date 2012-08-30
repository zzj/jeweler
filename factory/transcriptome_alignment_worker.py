import os
import json
import sys
import traceback
import argparse
from subprocess import call

def run_command(fd,comm):
	print(comm)
        ##os.system(comm)


def transcriptome_alignment_worker(args, refidtable, reffiletable):
    files=open(args.filelist).readlines()
    result_folder = 'result/'+os.path.basename(args.filelist)+'/transcriptome/'
    index_folder = result_folder+'/index/'
    cufflinks_folder = 'result/'+os.path.basename(args.filelist)+args.cufflinks_folder
    playpen_index='data/index/'

    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    ref_map={'F':'CAST','G':'PWK','H':'WSB'}
        
    for f in files:

        info = os.path.basename(f.strip().replace('.bam',''))
        resultsubfolder=result_folder+'/'+os.path.basename(f.strip().replace('.bam','')) + '/unmapped_aligned_to_' + ref_map[info[0]]
        if not os.path.exists(resultsubfolder):
            os.makedirs(resultsubfolder)
        unmapped_left_fastq_file = f.strip().replace('_merged.bam', '.left.unmapped.fastq')
        unmapped_right_fastq_file = f.strip().replace('_merged.bam', '.right.unmapped.fastq')
        gtf_input_file=	cufflinks_folder+'/'+os.path.basename(f.strip().replace('.bam',''))+'/transcripts.gtf'
        run_command(sys.stdout,'tophat  -o '+resultsubfolder+ '  -N 8 -p 4 -r 100 -T -G '+gtf_input_file+ ' --transcriptome-index '+ index_folder + ' '+
                    playpen_index+ref_map[info[0]]+'.fa ' + 
                    (unmapped_left_fastq_file)+' '+unmapped_right_fastq_file)
        if (info[0] == info[1]):
            resultsubfolder=result_folder+'/'+os.path.basename(f.strip().replace('.bam','')) + '/unmapped_aligned_to_' + ref_map[info[1]]
            if not os.path.exists(resultsubfolder):
                os.makedirs(resultsubfolder)
            run_command(sys.stdout, 'tophat  -o '+resultsubfolder+ ' -x 100 -N 8 -n 8 -p 4 -r 100 -T -G '+gtf_input_file+ ' --transcriptome-index '+ index_folder + ' '+
                        playpen_index+ref_map[info[1]]+'.fa ' + unmapped_left_fastq_file + ' ' + unmapped_right_fastq_file)
            