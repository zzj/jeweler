## get jeweler reuslt from the folder

find resfind result/merged_list/jeweler/ -maxdepth 2 -name "*merged.info"  | sort > info/merged_jeweler_result


## get all fasta files 

find . -name "*_sequence.txt" |sort > ~/Research/rna_seq/jeweler/info/fastq  


