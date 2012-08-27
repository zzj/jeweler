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


## TODO: merge left and right fastq file, find intersection.
## TODO: merge F1 fastq file, find intersection.
fastq_left = dict()
for infile in glob.glob(os.path.join(bamfolder + "*.left.fastq")):
        basename = os.path.basename(infile)
        index = basename.find('_aligned_to_')
        key = basename[:index]
        if (key not in fastq_left):
                fastq_left[key]=[]
        fastq_left[key].append(infile)

fastq_right = dict()
for infile in glob.glob(os.path.join(bamfolder + "*.right.fastq")):
        basename = os.path.basename(infile)
        index = basename.find('_aligned_to_')
        key = basename[:index]
        if (key not in fastq_right):
                fastq_right[key]=[]
        fastq_right[key].append(infile)

for key, value in table.items():
        if (len(value)==2):
                print('~/bin/bin/python2.7 mergeBam.py ' + ' '.join(value) + ' '+output_folder + '/'+key+'_merged.bam ')
                # if (key in fastq_left):
                #         print('mv '+ (fastq_left[key][0]) + ' ' + output_folder + '/' + key +".left.unmapped.fastq")
                #         print('mv '+ (fastq_right[key][0]) + ' ' + output_folder + '/' + key +".right.unmapped.fastq")
        else:
                continue
                bn = os.path.basename(value[0])
                if (bn[0] == bn[1]):
                        print('cp ' + value[0] + ' ' + output_folder +'/')
                else :
                        print('Error ' + value[0] )