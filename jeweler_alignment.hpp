#ifndef _JEWELER_ALIGNMENT_HPP_
#define _JEWELER_ALIGNMENT_HPP_
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <vector>
#include "proto/jeweler.pb.h"
#include "common.hpp"

using namespace std;
using namespace BamTools;

class SewingMachine;
class JewelerAlignment;

class AlignmentExpert {
public:
    virtual void initialize(JewelerAlignment *al);

    virtual void study_matched_seq(JewelerAlignment *al, int genome_start,
                                   int alignment_start, int length) ;

    virtual void study_only_read_seq(JewelerAlignment *al, int genome_start,
                                     int alignment_start, int length) ;

    virtual void study_only_genome_seq(JewelerAlignment *al, int genome_start,
                                       int alignment_start, int length) ;

    virtual void study_neither_exist_seq(JewelerAlignment *al, int genome_start,
                                          int alignment_start, int length);
};

class AlignmentExpertStarter : public AlignmentExpert {
public:

    virtual void initialize(JewelerAlignment *al);


    virtual void study_matched_seq(JewelerAlignment *al,
                           int genome_start,
                           int alignment_start,
                           int length);

    virtual void study_only_read_seq(JewelerAlignment *al,
                             int genome_start,
                             int alignment_start,
                             int length);
};

class JewelerAlignment : public BamAlignment {
public:
    _PROXY_(JewelerAlignment)
	// if two reads are merged, record the read position
	vector<int> genome_position;
    
    read_only<bool> is_multiple_alignment;

    void dump_data(Jeweler::EarringsData::Read *read);
    void set_is_multiple_alignment(const bool a);
    void set_is_multiple_alignment(const SewingMachine *sm);
    void jeweler_initialize(const SewingMachine *sm);
    void investigate(AlignmentExpert *ae);
    int GetStartPosition();
    int GetEndPosition();
    void glue(JewelerAlignment *that);
    int get_skipped_region(int skipped_alignment_length,
                           vector<CigarOp> &cigar_data,
                           int &skipped_length);

private:
    int is_first_ahead;
    int skip_length;
};

int get_read_position(JewelerAlignment *al, const int i);



#endif
