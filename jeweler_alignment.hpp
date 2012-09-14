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

class JewelerAlignment : public BamAlignment {
public:
    _PROXY_(JewelerAlignment)
	// if two reads are merged, record the read position
	vector<int> read_position;
	map<int, int> genome_position;
    void dump_data(Jeweler::EarringsData::Read *read);
    void set_is_multiple_alignment(const bool a);
    void set_is_multiple_alignment(const SewingMachine *sm);
    void jeweler_initialize(const SewingMachine *sm);
    read_only<bool> is_multiple_alignment;
private:

};

int get_read_position(JewelerAlignment *al, const int i);
#endif
