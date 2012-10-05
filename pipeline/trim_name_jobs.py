import glob
import os

for infile in glob.glob(os.path.join("../data/simulation_bam_merge_mapsplice/", "*_merged.bam")):
    print("python trim_name.py " + infile)