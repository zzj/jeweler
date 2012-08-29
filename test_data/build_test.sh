bowtie-build test_genome.fa test_genome
python generate_read.py
tophat test_genome left.fq right.fq 
