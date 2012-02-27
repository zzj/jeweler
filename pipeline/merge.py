import os, sys, time

import glob

bamfolder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new_combined/'
output_folder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new_merged/'

ref_map={'F':'CAST','G':'PWK','H':'WSB'}

table = dict()
os.system('mkdir '+output_folder)
for infile in glob.glob(os.path.join(bamfolder + "*.bam")):
        basename = os.path.basename(infile)
        index = basename.find('_aligned_to_')
        key = basename[:index]
        if (key not in table):
                table[key]=[]
        table[key].append(infile)

for key, value in table.items():
        if (len(value)==2):
                print('~/bin/bin/python2.7 mergeBam.py ' + ' '.join(value) + ' '+output_folder + '/'+key+'_merged.bam ')
        else:
                bn = os.path.basename(value[0])
                if (bn[0] == bn[1]):
                        print('cp ' + value[0] + ' ' + output_folder +'/')
                else :
                        print('Error ' + value[0] )