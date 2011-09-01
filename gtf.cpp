#include "gtf.hpp"
using namespace std;

transcript::transcript(){
	
}

int load_gtf_file(string gtf_filename, vector<transcript> &transcripts){


	FILE *fd=file_open(gtf_filename.c_str(),"r");
	char * temp=(char *)malloc(MAXLINE+1);
	char * line;
	char gene_id[100];
	char transcript_id[100];
	transcript t=transcript();
	int num=-1;


	transcripts.clear();
	while(fgets(temp, MAXLINE, fd)!=NULL){
		vector<string> strs;
		line=trim(temp);
		boost::split(strs,line,boost::is_any_of("\t"));
		if (strs.size()!=9){
			fprintf(stderr, "ERROR: GTF file %s corrupted \n",
					gtf_filename.c_str());
			exit(0);
		}
		if (strs[2]=="transcript" && num != -1){
			transcripts.push_back(t);
			t=transcript();
			num++;
		}
		else if (strs[2]=="transcript"){
			num++;
		}
		else if (strs[2]=="exon"){
			t.exon_start.push_back(atoi(strs[3].c_str()));
			t.exon_end.push_back(atoi(strs[4].c_str()));
		}
		sscanf(strs[8].c_str(),"gene_id \"%[^\"]\"; transcript_id \"%[^\"]\";",gene_id,transcript_id);

		t.name=transcript_id;

	}
	
	transcripts.push_back(t);

	free(temp);
	fclose(fd);
	return 0;
}
