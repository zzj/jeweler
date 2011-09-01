CC=-O2

all: jeweler

jeweler: jeweler.o common.o fasta.o gtf.o rna_read.o
	g++ $(CC) jeweler.o common.o fasta.o gtf.o rna_read.o -o jeweler

jeweler.o: jeweler.cpp 
	g++  $(CC) -c jeweler.cpp

common.o: common.cpp common.hpp
	g++  $(CC) -c common.cpp

fasta.o: fasta.cpp fasta.hpp
	g++ $(CC) -c fasta.cpp

gtf.o: gtf.cpp gtf.hpp
	g++ $(CC) -c gtf.cpp

rna_read.o: rna_read.cpp rna_read.hpp
	g++ $(CC) -c rna_read.cpp

clean: 
	rm *.o
	rm jeweler

test: all
	./jeweler -i test.info

deploy:
	cp scripts/* /playpen/rna_seq_comparison/
