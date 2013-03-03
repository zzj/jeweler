import os, sys, socket, time


def run_command(fd,comm):
    print(comm)

method = "mapsplice"
method = "contextmap"
method = "tophat"
method = "tophat_pseudo"
simulation_root = "/lustre/scr/z/z/zzj/RNAseqSim/output/"
fq_folder = "/real/"
if method == "tophat":
    output_root = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam/"
    output_root1 = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge/"
elif method == "tophat_pseudo":
    output_root = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_pseudo/"
    output_root1 = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge_pseudo/"
    fq_folder = "/pseudo/"
elif method =="mapsplice":
    output_root = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_mapsplice/"
    output_root1 = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge_mapsplice/"
elif method =="contextmap":
    output_root = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_contextmap/"
    output_root1 = "/lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge_contextmap/"
if not os.path.exists(output_root):
    os.makedirs(output_root)
if not os.path.exists(output_root1):
    os.makedirs(output_root1)

ref_map={'F':'CAST','G':'PWK','H':'WSB'}
n = 20
for i in range(0, n):
    for c, f in ref_map.items():
        sample_id = c+c+"_"+format(i, '04d')+"_"+"F"
        output_folder = output_root + "/"+sample_id
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        input_file1 = simulation_root + f + fq_folder + str(i) + "/" + "1.fq"
        input_file2 = simulation_root + f + fq_folder + str(i) + "/" + "2.fq"
        abundance_file = simulation_root + f + fq_folder + "abundance.txt." + str(i)
        if method == "tophat" or method == "tophat_pseudo":
            # print("tophat -o "+output_folder+"  -N 8 -p 4 -r 100 /lustre/scr/z/z/zzj/jeweler/data/index/"+f+".fa " + input_file1 + " " + input_file2)
            pass
        elif method == "mapsplice":
            # print("python ../../MapSplice_1.15.2/bin/mapsplice_segments.py --pairend -Q fq  -o "+output_folder + " -c /lustre/scr/z/z/zzj/jeweler/data/index/" + f  + "/  -B /lustre/scr/z/z/zzj/jeweler/data/index/"+f+".fa  -u " + input_file1 + "," + input_file2)
            # print("samtools view -bT ../data/index/" + f + ".fa " + output_folder + "/alignments.sam -o " + output_folder + "/accepted_hits_unsorted.bam ")
            # print("samtools sort -m 5000000000 "+ output_folder + "/accepted_hits_unsorted.bam " + output_folder + "/accepted_hits")
            pass
        elif method == 'contextmap':
            input_fa1 = simulation_root + f + fq_folder + str(i) + "/" + "1.fa"
            input_fa2 = simulation_root + f + fq_folder + str(i) + "/" + "2.fa"
            input_all = simulation_root + f + fq_folder + str(i) + "/" + "all.fa"
            # print("sed '/^@/!d;s//>/;N' " + input_file1 + " > " + input_fa1)
            # print("sed '/^@/!d;s//>/;N' " + input_file2 + " > " + input_fa2)
            # print("cat " + input_fa1 + " " + input_fa2 + " > " +input_all)
            # print("java -Xms7500m -Xmx15000m -XX:+UseConcMarkSweepGC -XX:NewSize=300M -XX:MaxNewSize=300M -XX:CMSInitiatingOccupancyFraction=40 -jar ContextMap_v1.1.5.jar mapper -reads "+ input_all + " -sam /nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/simulation_bam_mapsplice/" + sample_id + "/alignments.sam -readlen 100 -o " + output_folder + "/ -bwtbuild /nas02/apps/bowtie-0.12.8/bin/bowtie-build -bwtindex /nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/index/" + f + ".fa -genome /nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/index/" + f + "/ -rmap /nas02/home/z/z/zzj/Research/rna_seq/contextmap/rmap_x86_64")
            # print("\"samtools view -bT ../data/index/" + f + ".fa " + output_folder + "/mapping.sam -o " + output_folder + "/accepted_hits_unsorted.bam && samtools sort -m 5000000000 "+ output_folder + "/accepted_hits_unsorted.bam " + output_folder + "/accepted_hits\"")


        print("mv "+output_folder+"/accepted_hits.bam " + output_root1+"/"+sample_id+"_merged.bam")
        print("sort " + abundance_file + "> "+abundance_file+ ".sorted")
        print("join -2 2 "+abundance_file+".sorted" + " gid-tid-mapping.txt > "+output_root1+"/"+sample_id+".abundance")
