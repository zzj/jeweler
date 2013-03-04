import os, sys, time

import glob
import config

root = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/'
real_folder = root + 'simulation_bam' + config.tail + "_merge/"
pseudo_folder = root + 'simulation_bam' + config.tail + "_merge_pseudo/"
output_folder = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/simulation_bam_rp' + config.tail + '/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

simulation_root = "/lustre/scr/z/z/zzj/RNAseqSim/output/"
for infile in glob.glob(os.path.join(real_folder + "*merged.bam")):
    basename = os.path.basename(infile)

    # print('samtools merge ' + output_folder + basename + ' '
    #       + real_folder + basename + ' ' + pseudo_folder + basename + ' -f')
for infile in glob.glob(os.path.join(real_folder + "*.abundance")):
    basename = os.path.basename(infile)
    print("cat " + real_folder + basename + " " + pseudo_folder + basename +
          " > " + output_folder + basename)

    print("cp " + pseudo_folder + basename +
          "  " + output_folder + basename + ".pseudo")

    #print("sort " + abundance_file + "> "+abundance_file+ ".sorted")
    #print("join -2 2 "+abundance_file+".sorted" + " gid-tid-mapping.txt > "+output_root1+"/"+sample_id+".abundance")
