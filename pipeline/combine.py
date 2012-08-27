import os, sys, time

import glob

bamfolder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new/'
output_folder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new_combined/'

ref_map={'F':'CAST','G':'PWK','H':'WSB'}

table = dict()
left_unmapped = dict()
right_unmapped = dict()
pre_command = dict()
for infile in glob.glob(os.path.join(bamfolder + "*lane*")):
    basename = os.path.basename(infile)
    index = basename.find('_lane')
    key = basename[:index]
    if (key not in table):
        table[key]=[]
        pre_command[key]=[]
        left_unmapped[key] = []
        right_unmapped[key] = []

    table[key].append(infile+'/accepted_hits.bam')
    pre_command[key].append(('gunzip -c ' + infile+'/unmapped_left.fq.z > ' + infile+'/unmapped_left.fq') + " && " + ('gunzip -c ' + infile+'/unmapped_right.fq.z > '+ infile+'/unmapped_right.fq'))
    left_unmapped[key].append(infile+'/unmapped_left.fq')
    right_unmapped[key].append(infile+'/unmapped_right.fq')
        
for key, value in table.items():
    if len(value) != 4:
        continue
    print('samtools merge ' + output_folder + '/'+key+'.bam ' + ' '.join(value) + " -f")
        # left_join_command = "cat " + " ".join(left_unmapped[key]) + " > " + output_folder+'/'+key+".left.fastq"
        # right_join_command = "cat " + " ".join(right_unmapped[key]) + " > " + output_folder+'/'+key+".right.fastq"
        # fastq_command = " && ".join(pre_command[key]) + " && " + left_join_command + " && " + right_join_command
        # print(fastq_command)
        # if (len(value) != 4):
        #         print("Error:" + str(len(value)) + ' \nError: '.join(value), file = sys.stderr)
