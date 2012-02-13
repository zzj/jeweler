#include "gtf.hpp"
using namespace std;


int load_gtf_file(string gtf_filename, vector<Transcript *> &transcripts){


	FILE *fd=file_open(gtf_filename.c_str(),"r");
	char * temp=(char *)malloc(MAXLINE+1);
	char * line;
	char gene_id[100];
	char transcript_id[100];
	string chr;
	Transcript* t=new Transcript();
	int num=-1;

	int t_start = 0;
	int t_end = 0;

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
			t=new Transcript();
			// added by Weibo
			t_start = atoi(strs[3].c_str());
			t_end = atoi(strs[4].c_str());
			chr = strs[0];
			num++;
		}
		else if (strs[2]=="transcript"){
			num++;
			t_start = atoi(strs[3].c_str());
			t_end = atoi(strs[4].c_str());
			chr = strs[0];
		}
		else if (strs[2]=="exon"){
			int e_start = atoi(strs[3].c_str());
			int e_end = atoi(strs[4].c_str());
			t->exon_start.push_back(e_start);
			t->exon_end.push_back(e_end);
			// added by Weibo, add genome position
			for(int i = e_start; i <= e_end; ++i)
			{
				t->genome_pos.push_back(i);
			}
		}
		sscanf(strs[8].c_str(),"gene_id \"%[^\"]\"; transcript_id \"%[^\"]\";",gene_id,transcript_id);

		t->transcript_id=transcript_id;
		t->gene_id=gene_id;
		t->chr = chr;

	}

	// check error
	if (t->exon_start[0] != t_start)
	{
		fprintf(stderr, "%s:%d:%d\n",t->transcript_id.c_str(), t->exon_start[0],t_start);
		fprintf(stderr, "The start position of the first exon does not equal to the start position of the transcript \n");
		exit(0);
	}

	transcripts.push_back(t);

	free(temp);
	fclose(fd);
	return 0;
}
