#include "bam_info.hpp"
void BamInfo::initialize(BamReader &reader){
	this->sam_header = reader.GetHeader();
	this->references = reader.GetReferenceData();
}
