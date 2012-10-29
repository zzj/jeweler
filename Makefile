CC=g++
LIB=-Llib/bamtools/lib/ -L$(HOME)/bin/lib -Llib/leveldb-1.5.0/ -Ilib/glog-0.3.2/.lib/ 
LDFLAGS= -lz lib/bamtools/lib/libbamtools.a -lboost_filesystem -lprotobuf lib/leveldb-1.5.0/libleveldb.a lib/glog-0.3.2/.libs/libglog.a
CFLAGS= -g -O3 -std=gnu++0x -Wall -Wextra
INC=-Ilib/bamtools/include/ -Ilib/fastahack/ -I$(HOME)/bin/include -Ilib/leveldb-1.5.0/include -Ilib/glog-0.3.2/src/ -Ilib/protobuf-2.4.1/src

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

JEWELER_LIB_SOURCES =  jeweler_info.cpp jeweler_alignment.cpp read_matcher.cpp \
                 pileup.plot.cpp alignment_glue.cpp aligner.cpp transcript.cpp \
                 transcript_mismatcher.cpp earrings.cpp bracelet.cpp common.cpp \
                 fasta.cpp gtf_info.cpp gtf.cpp rna_read.cpp \
                 proto/jeweler.pb.cpp laboratory/cigar_holder.cpp \
                 laboratory/sewing_machine.cpp laboratory/bam_info.cpp \
                 graph/exon_node.cpp graph/graph.cpp  graph/path.cpp \
                 math/probability.cpp lib/fastahack/Fasta.cpp \
                 lib/fastahack/split.cpp

JEWELER_SOURCES= jeweler.cpp $(JEWELER_LIB_SOURCES)

TEST_JEWELER_SOURCES= test_common.cpp test_jeweler_alignment.cpp test_gtf.cpp \
                      test_jeweler_info.cpp test_transcript.cpp test_bracelet.cpp \
                      test_zleveldb.cpp laboratory/test_sewing_machine.cpp \
                      laboratory/test_cigar_holder.cpp test_alignment_glue.cpp

TEST_JEWELER_OBJECTS=$(TEST_JEWELER_SOURCES:.cpp=.o)
TEST_JEWELER = test_jeweler

JEWELER_EXECUTABLE=jeweler
JEWELER_OBJECTS= $(JEWELER_SOURCES:.cpp=.o)
JEWELER_LIB_OBJECTS=$(JEWELER_LIB_SOURCES:.cpp=.o)

APPRAISER_SOURCES= $(APPRAISER_CPP_SOURCES) $(APPRAISER_C_SOURCES)
APPRAISER_CPP_SOURCES=laboratory/appraiser.cpp  common.cpp  laboratory/metabam.cpp \
                      laboratory/bam_info.cpp  laboratory/sewing_machine.cpp \
                      proto/jeweler.pb.cpp \
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
OBJECTS=$(JEWELER_OBJECTS) $(APPRAISER_OBJECTS) $(TEST_BAM_OBJECTS) $(TEST_JEWELER_OBJECTS)

TESTS = $(TEST_JEWELER)

all: $(SOURCES) $(EXECUTABLE) $(TESTS)

$(JEWELER_EXECUTABLE): $(JEWELER_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

$(TEST_JEWELER): $(JEWELER_LIB_OBJECTS) $(TEST_JEWELER_OBJECTS) gtest_main.a 
	$(CC) $(LIB) $^ $(LDFLAGS) -lpthread  -o $@
	./$(TEST_JEWELER)

$(APPRAISER_EXECUTABLE): $(APPRAISER_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

$(TEST_BAM_EXECUTABLE): $(TEST_BAM_OBJECTS)
	$(CC) $(LIB) $^ $(LDFLAGS) -o $@

proto/jeweler.pb.cc: proto/jeweler.proto
	cd proto &&	protoc --cpp_out=. --python_out=../shop/ jeweler.proto

proto/jeweler.pb.cpp: proto/jeweler.pb.cc proto/jeweler.pb.h
	cp proto/jeweler.pb.cc proto/jeweler.pb.cpp

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

