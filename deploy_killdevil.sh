mkdir /lustre/scr/z/z/zzj/jeweler
mkdir /lustre/scr/z/z/zzj/jeweler/data
mkdir /lustre/scr/z/z/zzj/jeweler/data/index
mkdir /lustre/scr/z/z/zzj/jeweler/data/cegs_rnaseq_inbreds
mkdir /lustre/scr/z/z/zzj/jeweler/data/cegs_rnaseq_merged
mkdir /lustre/scr/z/z/zzj/jeweler/data/database/
ln -s /lustre/scr/z/z/zzj/jeweler/data data


scp -r zzj@csbio-desktop107.cs.unc.edu:/home/zzj/Research/rna_seq/data/ensembl/\* data/database/
scp -r zzj@csbio-desktop107.cs.unc.edu:/home/zzj/Research/rna_seq/data/genomes/\* data/index/





