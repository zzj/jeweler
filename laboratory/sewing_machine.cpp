
#include "sewing_machine.hpp"

locator::locator(int ref_id, int position, string cigar_string){
	this->ref_id=ref_id;
	this->position=position;
	this->cigar_string=cigar_string;
}

SewingMachine::SewingMachine(){
}

SewingMachine::SewingMachine(BamReader &reader){
	BamInfo::initialize(reader);
}


int SewingMachine::add_alignment(BamAlignment &al){
	locator l(al.RefID,al.Position, get_cigar_string(al));
	string is_first;
	if (al.IsFirstMate()){
		is_first="\tT";
	}
	else {
		is_first="\tF";
	}
	seqs[al.Name+is_first].push_back(l);
	return 0;
}

int SewingMachine::output_alignment_map(FILE * file){

	auto i=seqs.begin(); 
	int j;
	for ( ; i!=seqs.end();i++){
		fprintf(file,"%s\t%zd\t",i->first.c_str(),i->second.size());
		for (j=0;j<i->second.size();j++){
			fprintf(file,"%s\t%d\t%s\t",
					references[i->second[j].ref_id].RefName.c_str(),
					i->second[j].position,
					i->second[j].cigar_string.c_str());
		}
		fprintf(file,"\n");
	}
}

int SewingMachine::output_multiple_alignment_map(FILE * file){

	auto i=seqs.begin(); 
	int j;
	for ( ; i!=seqs.end();i++){
		if (i->second.size()<2) continue;
		fprintf(file,"%s\t%zd\t",i->first.c_str(),i->second.size());
		for (j=0;j<i->second.size();j++){
			fprintf(file,"%s\t%d\t%s\t",
					references[i->second[j].ref_id].RefName.c_str(),
					i->second[j].position,
					i->second[j].cigar_string.c_str());
		}
		fprintf(file,"\n");
	}
}
