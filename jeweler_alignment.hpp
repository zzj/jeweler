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
    virtual void study_matched_seq(JewelerAlignment *al, int genome_start,
                                   int alignment_start, int length) ;

    virtual void study_only_read_seq(JewelerAlignment *al, int genome_start,
                                     int alignment_start, int length) ;

    virtual void study_only_genome_seq(JewelerAlignment *al, int genome_start,
                                       int alignment_start, int length) ;

    virtual void study_neither_exist_seq(JewelerAlignment *al, int genome_start,
                                          int alignment_start, int length);
};


class JewelerAlignment : public BamAlignment {
public:
    _PROXY_(JewelerAlignment)
	// if two reads are merged, record the read position
	vector<int> read_position;
	map<int, int> genome_position;
    read_only<bool> is_multiple_alignment;

    void dump_data(Jeweler::EarringsData::Read *read);
    void set_is_multiple_alignment(const bool a);
    void set_is_multiple_alignment(const SewingMachine *sm);
    void jeweler_initialize(const SewingMachine *sm);
    void investigate(AlignmentExpert *ae);
    int GenomeStartPosition();
private:

};

int get_read_position(JewelerAlignment *al, const int i);



#endif
