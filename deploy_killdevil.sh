mkdir /lustre/scr/z/z/zzj/jeweler
mkdir /lustre/scr/z/z/zzj/jeweler/data
mkdir /lustre/scr/z/z/zzj/jeweler/data/index
mkdir /lustre/scr/z/z/zzj/jeweler/data/cegs_rnaseq_inbreds
mkdir /lustre/scr/z/z/zzj/jeweler/data/cegs_rnaseq_merged
mkdir /lustre/scr/z/z/zzj/jeweler/data/database/
ln -s /lustre/scr/z/z/zzj/jeweler/data data


scp -r zzj@csbio-desktop107.cs.unc.edu:/home/zzj/Research/rna_seq/jeweler/data/index/* data/index
scp -r zzj@csbio-storage001.cs.unc.edu:/csbiodataxw/RNAseq-nobackup4/zzj/jeweler/data/ data
##scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodataxu/RNAseq-nobackup3/cegs_rnaseq_merged  data
scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodataxu/RNAseq-nobackup3/cegs_rnaseq_inbreds  data
scp -r zzj@csbio-desktop107.cs.unc.edu:/home/zzj/Research/rna_seq/data/ensembl/* data/database/
scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodata/RNAseq-nobackup4/zzj/jeweler/data/index/* data/index/
scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodataxz/RNAseq-nobackup2/cegs_rnaseq_output2/combined/111031_111104/ data/cegs_rnaseq_merged/
mv data/cegs_rnaseq_merged/FF* data/cegs_rnaseq_inbreds/
mv data/cegs_rnaseq_merged/GG* data/cegs_rnaseq_inbreds/
mv data/cegs_rnaseq_merged/HH* data/cegs_rnaseq_inbreds/


##copy fastq files

scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodata/RNAseq-nobackup/cegs_rnaseq_111030/ data/fastq/
scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodata/RNAseq-nobackup/cegs_rnaseq_111104/ data/fastq/

scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodata/RNAseq-nobackup4/CEGS_3x3/RNAseq/1.fastq/111108_UNC10 data/fastq/
scp -r zzj@csbio-desktop107.cs.unc.edu:/csbiodata/RNAseq-nobackup4/CEGS_3x3/RNAseq/1.fastq/111123_UNC9 data/fastq/





