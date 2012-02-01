#ifndef _JEWELER_ALIGNMENT_HPP_
#define _JEWELER_ALIGNMENT_HPP_
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <vector>
using namespace std;
using namespace BamTools;
class JewelerAlignment : public BamAlignment {
public:
	// if two reads are merged, record the read position
	vector<int> read_position;
};
#endif
