#ifndef _JEWELER_ALIGNMENT_HPP_
#define _JEWELER_ALIGNMENT_HPP_
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <vector>
#include "proto/jeweler.pb.h"
#include "laboratory/cigar_holder.hpp"

using namespace std;
using namespace BamTools;
class JewelerAlignment : public BamAlignment {
public:
	// if two reads are merged, record the read position
	vector<int> read_position;
	vector<int> genome_position;
    void dump_data(Jeweler::EarringsData::Read *read);
    void set_is_multiple_alignment(bool a);
private:
    bool is_multiple_alignment;

};

int get_read_position(JewelerAlignment *al, int i);
#endif
