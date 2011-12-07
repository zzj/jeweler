
#include "transcript_info.hpp"

transcript_info::transcript_info(string gene_id,string folder,string gtf_filename,
								 string paternal_seq_filename,
								 string maternal_seq_filename,
								 string bam_read_filename){
	this->gene_id=gene_id;
	this->folder=folder;
	this->gtf_filename=gtf_filename;
	this->paternal_seq_filename=paternal_seq_filename;
	this->maternal_seq_filename=maternal_seq_filename;
	this->read_seq_filename=bam_read_filename;
}

