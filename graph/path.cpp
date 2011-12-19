#include "path.hpp"

Path::Path(vector<ExonNode*> &current_path){
	path=current_path;
}


bool Path::is_mirrored(Path &a){
	// must be informative path (at least contain one allele specific
	// exon)
	
	bool ret=false;
	if (a.path.size()!=path.size()){
		return false;
	}
	
	for (int i=0;i<path.size();i++){
		if (!path[i]->is_mirrored(a.path[i])){
			return false;
		}
		if (path[i]->origin!=EXON_NO_INFO){
			ret=true;
		}
	}
	return ret;
}

bool Path::is_informative(){
	for (int i=0;i<path.size();i++){
		if (path[i]->origin!=EXON_NO_INFO){
			return true;
		}
	}
	return false;
}

bool Path::is_equal(Path &a){
	for (int i=0;i<path.size();i++){
		if (!path[i]->is_equal(a.path[i])){
			return false;
		}
	}
	return true;
}


int Path::dump_path(FILE * file){
	int i;
	if (is_valid()){
		fprintf(file,"%s ",path[0]->detach().c_str());
		for (i=1;i<path.size();i++){
			fprintf(file,"-> %s ",path[i]->detach().c_str());
		}
		fprintf(file,"\n");
	}
}

bool Path::is_valid(){
	for (int i=0;i<path.size();i++){
		if (path[i]==NULL){
			return false;
		}
	}
	return true;
}

bool Path::is_compatible(){
	bool maternal=false, paternal=false;
	
	for (int i=0;i<path.size();i++){
		if (path[i]==NULL){
			continue;
		}
		if (path[i]->origin==EXON_PATERNAL){
			paternal=true;
		}
		if (path[i]->origin==EXON_MATERNAL){
			maternal=true;
		}
	}
	return !(paternal && maternal);
	
}
