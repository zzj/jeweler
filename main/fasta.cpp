#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include "fasta.hpp"
#include "common.hpp"
using namespace std;



seq_read::seq_read(string name, string seq) {
	this->name=name;
	this->seq=seq;
	size=this->seq.size();
}
bool operator < (const seq_read &a, const seq_read &b) {
	return a.name<b.name;

}

int load_fasta_file(string fasta_filename, vector<seq_read *> &seq_reads) {
	FILE *fd=file_open(fasta_filename.c_str(),"r");
	string name;
	string seq;
	char * temp=(char *)malloc(MAXLINE+1);
	char * line;

	seq_reads.clear();
	
	while(fgets(temp, MAXLINE, fd)!=NULL) {
		line=trim(temp);
		if (strlen(line)==0 || line[0]=='>') {
			if (seq.size()>0) {
				seq_read * r=new seq_read(name,seq);
				seq_reads.push_back(r);
			}
			seq.clear();
			strtok(line," \t\n");
			name=line+1;
		}
		else {
			seq=seq+line;
		}
	}
	if (seq.size()>0) {
		seq_read * r=new seq_read(name,seq);
		seq_reads.push_back(r);
	}
	fclose(fd);
	free(temp);
	return 0;
}

void write_fasta_file(const string fasta_filename, const string &name, const string &seq) {
	FILE *fd = file_open(fasta_filename.c_str(), "w+");
	write_fasta_file(fd, name,seq);
	fclose(fd);
}
void write_fasta_file(FILE *fd, const string &name, const string &seq) {
	fprintf(fd,">%s\n", name.c_str());
	for (size_t i = 0; i < seq.size(); i = i + 80) {
		fprintf(fd,"%s\n", seq.substr(i, 80).c_str());
	}
}
