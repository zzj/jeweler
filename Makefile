all: jeweler


jeweler: jeweler.o common.o fasta.o gtf.o rna_read.o
	g++ jeweler.o common.o fasta.o gtf.o rna_read.o -o jeweler

jeweler.o: jeweler.cpp 
	g++ -c jeweler.cpp

common.o: common.cpp common.hpp
	g++ -c common.cpp

fasta.o: fasta.cpp fasta.hpp
	g++ -c fasta.cpp

gtf.o: gtf.cpp gtf.hpp
	g++ -c gtf.cpp

rna_read.o: rna_read.cpp rna_read.hpp
	g++ -c rna_read.cpp

clean: 
	rm *.o
	rm jeweler

test: all
	./jeweler -i test.info

deploy:
	cp scripts/* /playpen/rna_seq_comparison/
