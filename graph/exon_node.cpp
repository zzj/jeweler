#include "exon_node.hpp"

ExonNode::ExonNode(int start, int end, int origin){
	this->start=start;
	this->end=end;
	this->origin=origin;
}

string ExonNode::detach(){
	string ret;
	char temp[1000];
	sprintf(temp,"%d:%d:%d",start, end, origin);
	ret=temp;
	return ret;
}
