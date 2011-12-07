CC=-g -std=gnu++0x
INC=-Ilib/bamtools/include/ -Ilib/fastahack/
LIB=-Llib/bamtools/lib/


all: jeweler appraiser

jeweler: jeweler.o earrings.o common.o fasta.o gtf.o rna_read.o transcript_info.o transcript.o cigar_holder.o landscape.plot.o
	g++ $(CC) jeweler.o earrings.o common.o fasta.o gtf.o rna_read.o transcript_info.o  landscape.plot.o transcript.o cigar_holder.o -lz lib/bamtools/lib/libbamtools.a -o jeweler



jeweler.o: jeweler.cpp  jeweler.hpp
	g++  $(INC)  $(CC) -c jeweler.cpp 

transcript_info.o: transcript_info.cpp transcript_info.hpp
	g++ $(INC) $(CC) -c transcript_info.cpp

landscape.plot.o: landscape.plot.cpp landscape.plot.hpp
	g++ $(INC) $(CC) -c landscape.plot.cpp

transcript.o: transcript.cpp transcript.hpp
	g++ $(INC) $(CC) -c transcript.cpp

earrings.o: earrings.cpp earrings.hpp 
	g++  $(INC) $(CC) -c earrings.cpp	

common.o: common.cpp common.hpp
	g++  $(CC) -c common.cpp

fasta.o: fasta.cpp fasta.hpp
	g++ $(CC) -c fasta.cpp

gtf.o: gtf.cpp gtf.hpp
	g++ $(INC) $(CC) -c gtf.cpp

rna_read.o: rna_read.cpp rna_read.hpp
	g++ $(CC) -c rna_read.cpp

clean: 
	rm *.o
	rm jeweler appraiser

test: all
	./jeweler -i test.info

appraiser: appraiser.o common.o metabam.o bam_info.o sewing_machine.o cigar_holder.o Fasta.o split.o disorder.o 
	g++ $(INC) $(LIB) $(CC) appraiser.o sewing_machine.o metabam.o  bam_info.o cigar_holder.o  Fasta.o split.o disorder.o -lz lib/bamtools/lib/libbamtools.a common.o  -o appraiser 

test_bam: test_bam.cpp Fasta.o split.o
	g++  test_bam.cpp $(INC) $(LIB) $(CC) Fasta.o split.o cigar_holder.o -lz lib/bamtools/lib/libbamtools.a   -o test_bam

appraiser.o: laboratory/appraiser.cpp laboratory/appraiser.hpp
	g++ $(INC) $(CC) -c laboratory/appraiser.cpp  


metabam.o: laboratory/metabam.cpp laboratory/metabam.hpp
	g++ $(INC) $(CC) -c laboratory/metabam.cpp  

bam_info.o: laboratory/bam_info.cpp laboratory/bam_info.hpp
	g++ $(INC) $(CC) -c laboratory/bam_info.cpp  		

cigar_holder.o: laboratory/cigar_holder.cpp laboratory/cigar_holder.hpp
	g++ $(INC) $(CC) -c laboratory/cigar_holder.cpp  		


sewing_machine.o: laboratory/sewing_machine.cpp laboratory/sewing_machine.hpp
	g++ $(INC) $(CC) -c laboratory/sewing_machine.cpp  		

Fasta.o: lib/fastahack/Fasta.h lib/fastahack/Fasta.cpp
	g++ $(INC) $(CC) -c lib/fastahack/Fasta.cpp

split.o: lib/fastahack/split.h lib/fastahack/split.cpp
	g++ -c lib/fastahack/split.cpp

disorder.o: lib/fastahack/disorder.c lib/fastahack/disorder.h
	g++ -c lib/fastahack/disorder.c
