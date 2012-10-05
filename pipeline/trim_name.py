import pysam
import glob
import os
infile = os.sys.argv[1]
pysam.index(infile)
samfile = pysam.Samfile(infile, "rb")
new_samfile = pysam.Samfile(infile[:-4] + "_trim.bam", "wb", template=samfile)

for read in samfile.fetch():
    read.qname = read.qname[:-2]
    new_samfile.write(read)

samfile.close()
new_samfile.close()
