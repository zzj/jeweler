#include "bam_info.hpp"
int BamInfo::initialize(BamReader &reader){
	this->sam_header = reader.GetHeader();
	this->references = reader.GetReferenceData();
	return 0;
}
