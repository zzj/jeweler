CC=-g -std=gnu++0x
INC=-Ilib/bamtools/include/
LIB=-Llib/bamtools/lib/


all: jeweler appraiser

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
	rm jeweler appraiser

test: all
	./jeweler -i test.info

appraiser: appraiser.o common.o metabam.o bam_info.o sewing_machine.o cigar_holder.o
	g++ $(INC) $(LIB) $(CC) appraiser.o sewing_machine.o metabam.o  bam_info.o cigar_holder.o -lz lib/bamtools/lib/libbamtools.a common.o  -o appraiser 

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
