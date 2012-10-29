import glob
import os

for infile in glob.glob(os.path.join("../data/simulation_bam_merge_contextmap/", "*_merged.bam")):
    print("python trim_name.py " + infile)