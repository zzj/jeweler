#include "exon_node.hpp"
#include "main/constants.hpp"
#include "main/common.hpp"
#include "main/jeweler_alignment.hpp"

ExonNode::ExonNode(int start, int end, int origin){
	this->start=start;
	this->end=end;
	this->origin=origin;
}

void ExonNode::insert_reads( set<JewelerAlignment *> & reads){
	this->reads.insert(reads.begin(),reads.end());
}


void ExonNode::insert_allele_reads( set<JewelerAlignment *> & reads){
	this->allele_reads.insert(reads.begin(),reads.end());
}

string ExonNode::detach(){
	string ret;
	char temp[1000];
	sprintf(temp,"%d:%d:%zu:%zu:%d",start, end, reads.size(),allele_reads.size(), origin);
	ret=temp;
	return ret;
}

bool ExonNode::is_mirrored(ExonNode * a){
	if (a->start==this->start && a->end==this->end){
		if (a->origin==EXON_NO_INFO && this->origin==a->origin){
			return true;
		}
		else{
			if ((a->origin==EXON_PATERNAL && this->origin==EXON_MATERNAL) ||
				(a->origin==EXON_MATERNAL && this->origin==EXON_PATERNAL)){
				return true;
			}
			else return false;
		}
	}
	else return false;
}

bool ExonNode::is_equal(ExonNode * a){
	if (a->start==this->start && a->end==this->end && a->origin==this->origin){
		return true;
	}	
	return false;
}
