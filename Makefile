CFLAGS= -g -std=gnu++0x -Wall -Wextra
INC=-Ilib/bamtools/include/ -Ilib/fastahack/ -I$(HOME)/bin/include

# UNIT TEST MAKEFILE

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = lib/gtest-1.6.0/

# Where to find user code.
USER_DIR = .

# Flags passed to the preprocessor.
CXXFLAGS= $(CFLAGS)

# Flags passed to the preprocessor.
CPPFLAGS=$(INC) -I$(GTEST_DIR)/include

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
CC=g++
LIB=-Llib/bamtools/lib/ -L$(HOME)/bin/lib
LDFLAGS= -lz lib/bamtools/lib/libbamtools.a -lboost_filesystem

JEWELER_LIB_SOURCES =  jeweler_alignment.cpp read_matcher.cpp \
                 pileup.plot.cpp alignment_glue.cpp aligner.cpp transcript.cpp \
                 transcript_mismatcher.cpp earrings.cpp bracelet.cpp common.cpp \
                 fasta.cpp gtf.cpp rna_read.cpp laboratory/cigar_holder.cpp \
                 laboratory/sewing_machine.cpp laboratory/bam_info.cpp \
                 graph/exon_node.cpp graph/graph.cpp  graph/path.cpp \
                 math/probability.cpp lib/fastahack/Fasta.cpp \
                 lib/fastahack/split.cpp

JEWELER_SOURCES= jeweler.cpp $(JEWELER_LIB_SOURCES)

TEST_JEWELER_SOURCES= test_jeweler.cpp
TEST_JEWELER_OBJECTS=$(TEST_JEWELER_SOURCES:.cpp=.o)
TEST_JEWELER = test_jeweler

JEWELER_EXECUTABLE=jeweler
JEWELER_OBJECTS=$(JEWELER_SOURCES:.cpp=.o)
JEWELER_LIB_OBJECTS=$(JEWELER_LIB_SOURCES:.cpp=.o)

APPRAISER_SOURCES= $(APPRAISER_CPP_SOURCES) $(APPRAISER_C_SOURCES)
APPRAISER_CPP_SOURCES=laboratory/appraiser.cpp  common.cpp  laboratory/metabam.cpp \
                      laboratory/bam_info.cpp  laboratory/sewing_machine.cpp \
                      laboratory/cigar_holder.cpp  lib/fastahack/Fasta.cpp \
                      lib/fastahack/split.cpp

APPRAISER_C_SOURCES=lib/fastahack/disorder.c
APPRAISER_EXECUTABLE= appraiser
APPRAISER_CPP_OBJECTS=$(APPRAISER_CPP_SOURCES:.cpp=.o)
APPRAISER_C_OBJECTS=$(APPRAISER_C_SOURCES:.c=.o)
APPRAISER_OBJECTS=$(APPRAISER_C_OBJECTS) $(APPRAISER_CPP_OBJECTS)

TEST_BAM_SOURCES= test_bam.cpp lib/fastahack/Fasta.cpp \
                  lib/fastahack/split.cpp  laboratory/cigar_holder.cpp
TEST_BAM_EXECUTABLE=test_bam
TEST_BAM_OBJECTS=$(TEST_BAM_SOURCES:.cpp=.o)

EXECUTABLE=$(JEWELER_EXECUTABLE) $(APPRAISER_EXECUTABLE) $(TEST_BAM_EXECUTABLE)
SOURCES=$(JEWELER_SOURCES) $(APPRAISER_SOURCES)
OBJECTS=$(JEWELER_OBJECTS) $(APPRAISER_OBJECTS) $(TEST_BAM_OBJECTS)

all: $(SOURCES) $(EXECUTABLE) $(TEST)

$(JEWELER_EXECUTABLE): $(JEWELER_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

$(TEST_JEWELER): $(JEWELER_LIB_OBJECTS) $(TEST_JEWELER_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -lpthread  -o $@

$(APPRAISER_EXECUTABLE): $(APPRAISER_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

$(TEST_BAM_EXECUTABLE): $(TEST_BAM_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -MD -c -o $@ $<
	cp $*.d $*.P; \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
			-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
			rm -f $*.d

-include $(OBJECTS:.o=.P)


gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^



clean:
	rm $(sort $(OBJECTS)) $(EXECUTABLE) $(TESTS) gtest.a gtest_main.a

