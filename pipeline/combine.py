import os, sys, time

import glob

bamfolder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new/'
output_folder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new_combined/'

ref_map={'F':'CAST','G':'PWK','H':'WSB'}

table = dict()

for infile in glob.glob(os.path.join(bamfolder + "*lane*")):
        basename = os.path.basename(infile)
        index = basename.find('_lane')
        key = basename[:index]
        if (key not in table):
                table[key]=[]

        table[key].append(infile+'/accepted_hits.bam')



for key, value in table.items():
        print('samtools merge ' + output_folder + '/'+key+'.bam ' + ' '.join(value) + " -f")