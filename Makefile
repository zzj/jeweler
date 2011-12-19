CC=g++
CFLAGS=-g -std=gnu++0x
INC=-Ilib/bamtools/include/ -Ilib/fastahack/ -I$(HOME)/bin/include
LIB=-Llib/bamtools/lib/ -I$(HOME)/bin/lib
LDFLAGS= -lz lib/bamtools/lib/libbamtools.a 


JEWELER_SOURCES= jeweler.cpp transcript_info.cpp landscape.plot.cpp transcript.cpp earrings.cpp common.cpp fasta.cpp gtf.cpp rna_read.cpp laboratory/cigar_holder.cpp graph/exon_node.cpp graph/graph.cpp  graph/path.cpp
JEWELER_EXECUTABLE=jeweler 
JEWELER_OBJECTS=$(JEWELER_SOURCES:.cpp=.o) 

APPRAISER_SOURCES= $(APPRAISER_CPP_SOURCES) $(APPRAISER_C_SOURCES)
APPRAISER_CPP_SOURCES=laboratory/appraiser.cpp  common.cpp  laboratory/metabam.cpp  laboratory/bam_info.cpp  laboratory/sewing_machine.cpp  laboratory/cigar_holder.cpp  lib/fastahack/Fasta.cpp  lib/fastahack/split.cpp  
APPRAISER_C_SOURCES=lib/fastahack/disorder.c
APPRAISER_EXECUTABLE= appraiser
APPRAISER_CPP_OBJECTS=$(APPRAISER_CPP_SOURCES:.cpp=.o) 
APPRAISER_C_OBJECTS=$(APPRAISER_C_SOURCES:.c=.o) 
APPRAISER_OBJECTS=$(APPRAISER_C_OBJECTS) $(APPRAISER_CPP_OBJECTS)

TEST_BAM_SOURCES= test_bam.cpp lib/fastahack/Fasta.cpp  lib/fastahack/split.cpp  laboratory/cigar_holder.cpp
TEST_BAM_EXECUTABLE=test_bam
TEST_BAM_OBJECTS=$(TEST_BAM_SOURCES:.cpp=.o)

EXECUTABLE=$(JEWELER_EXECUTABLE) $(APPRAISER_EXECUTABLE) $(TEST_BAM_EXECUTABLE)
SOURCES=$(JEWELER_SOURCES) $(APPRAISER_SOURCES)
OBJECTS=$(JEWELER_OBJECTS) $(APPRAISER_OBJECTS) $(TEST_BAM_OBJECTS)


all: $(SOURCES) $(EXECUTABLE)

$(JEWELER_EXECUTABLE): $(JEWELER_OBJECTS) 
	$(CC) $(LIB) $(JEWELER_OBJECTS) $(LDFLAGS) -o $@

$(APPRAISER_EXECUTABLE): $(APPRAISER_OBJECTS) 
	$(CC) $(LIB) $(APPRAISER_OBJECTS) $(LDFLAGS) -o $@

$(TEST_BAM_EXECUTABLE): $(TEST_BAM_OBJECTS) 
	$(CC) $(LIB) $(TEST_BAM_OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(INC) $(CFLAGS) -MD -c -o $@ $<
	cp $*.d $*.P; \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
			-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
			rm -f $*.d

-include $(OBJECTS:.o=.P)


clean:
	rm $(sort $(OBJECTS)) $(EXECUTABLE)
	rm *.o

