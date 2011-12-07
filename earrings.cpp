#include "earrings.hpp"

Earrings::Earrings(transcript_info *info){
	load_transcript_data(info);	
}

int Earrings::load_transcript_data(transcript_info * info){
	int i,j;

	transcript_helper(info->paternal_seq_filename,info->gtf_filename, 
					  paternal_transcripts );
	transcript_helper(info->maternal_seq_filename,info->gtf_filename, 
					  maternal_transcripts );
	if (paternal_transcripts.size()!=maternal_transcripts.size() 
		|| paternal_transcripts.size()==0){
		fprintf(stderr, 
				"ERROR: number of transcripts does not match or no reads at all for gene %s  at %s:%d\n",
				info->gene_id.c_str(),__FILE__, __LINE__);
		exit(0);
	}

	for (i=0;i<paternal_transcripts.size();i++){
		transcript *p=paternal_transcripts[i];
		transcript *m=maternal_transcripts[i];
		if (p->seq.size()!=m->seq.size()){

			fprintf(stderr, 
					"ERROR: transcript sequence size does not match at gene %s  at %s:%d\n",
					p->name.c_str(),__FILE__, __LINE__);
			exit(0);
		}
		int num=0;

		// find SNP position by given both paternal and maternal transcripts sequences
		for (j=0;j<p->seq.size();j++){
			if (p->seq[j]!=m->seq[j]){
				p->snp_pos.push_back(j);
				m->snp_pos.push_back(j);
				p->alleles.push_back(p->seq[j]);
				m->alleles.push_back(m->seq[j]);
				num++;
			}
		}

	}
	return 0;
}

int Earrings::transcript_helper(string seq_filename,string gtf_file, 
							   vector<transcript *> &transcripts){
	//assuming gtf file has the same order of transcripts with the seq files
	int i,j;
	vector<seq_read*> sr;
	load_fasta_file(seq_filename,sr);
	load_gtf_file(gtf_file,transcripts);
	if (transcripts.size()!=sr.size()){
		fprintf(stderr, "ERROR: sequence file (%d sequences) does not match gtf file (%d transcripts) at %s:%d\n",
				(int)sr.size(),(int)transcripts.size(),__FILE__, __LINE__);
		exit(0);

	}

	// TODO: put these code into gtf.cpp files.
	// All transcript class operations should be done in transcripts. 
	for (i=0;i<sr.size();i++){
		for (j=0;j<transcripts.size();j++){
			if (transcripts[j]->name != sr[i]->name){
				continue;
			}
			transcripts[j]->seq=sr[i]->seq;
			transcripts[j]->noninformative_mismatches.resize(transcripts[j]->seq.size(),0);
			//check errors
			if (transcripts[j]->seq.size() != transcripts[j]->genome_pos.size())
			{
				fprintf(stderr, "something wrong in inferring the genome position");
				exit(0);
			}
			else
			{
				//fprintf(stderr, "size the transcirpt is %d bases\n", transcripts[j]->genome_pos.size());
				transcripts[j]->Anoninformative_mismatches.resize(transcripts[j]->genome_pos.size());
				transcripts[j]->Cnoninformative_mismatches.resize(transcripts[j]->genome_pos.size());
				transcripts[j]->Gnoninformative_mismatches.resize(transcripts[j]->genome_pos.size());
				transcripts[j]->Tnoninformative_mismatches.resize(transcripts[j]->genome_pos.size());
			}
		}
	}
	return 0;
}
