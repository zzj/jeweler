#include "earrings.hpp"

Earrings::Earrings(TranscriptInfo *info){

	this->info=info;
	// load maternal and paternal transcripts sequences.
	load_transcript_data(info);	
	// load bamalignment by bamtools.
	load_read_data(info);
	// align reads back to the transcripts
	align_reads();
	// build graph;
	build_graph();
	test_allele_specific_transcript();
}


Earrings::~Earrings(){
	int i=0;
	for (i=0;i<maternal_transcripts.size();i++){
		delete maternal_transcripts[i];
		delete paternal_transcripts[i];
	}
	for (i=0;i<bam_reads.size();i++){
		delete bam_reads[i];
	}
}

int Earrings::load_read_data(TranscriptInfo *info){
	string bam_filename=info->read_seq_filename;

	BamReader reader;

	if (!reader.Open(bam_filename)){
		fprintf(stderr,"Cannot open bam file!\n");
		exit(0);
	}
	BamAlignment *al=new BamAlignment();

	while(reader.GetNextAlignment(*al)){
		bam_reads.push_back(al);
		al=new BamAlignment();
	}
	delete al; // delete the last unused one
	
	return 0;
}



int Earrings::load_transcript_data(TranscriptInfo * info){
	int i,j,k;

	transcript_helper(info->maternal_seq_filename,info->gtf_filename, 
					  maternal_transcripts );
	transcript_helper(info->paternal_seq_filename,info->gtf_filename, 
					  paternal_transcripts );
	if (paternal_transcripts.size()!=maternal_transcripts.size() 
		|| paternal_transcripts.size()==0){
		fprintf(stderr, 
				"ERROR: number of transcripts does not match or no reads at all for gene %s  at %s:%d\n",
				info->gene_id.c_str(),__FILE__, __LINE__);
		exit(0);
	}

	for (i=0;i<paternal_transcripts.size();i++){
		Transcript *p=paternal_transcripts[i];
		Transcript *m=maternal_transcripts[i];
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
		p->allele_exon.resize(m->snp_pos.size());
		m->allele_exon.resize(m->snp_pos.size());
		p->num_info_reads_per_exon.resize(p->exon_start.size(),0);
		m->num_info_reads_per_exon.resize(m->exon_start.size(),0);
		p->origin=EXON_PATERNAL;
		m->origin=EXON_MATERNAL;
		p->num_alleles_per_exon.resize(p->exon_start.size(),0);
		m->num_alleles_per_exon.resize(p->exon_start.size(),0);

		for (j=0;j<m->snp_pos.size();j++){
			// find the exon
			int exon_id=-1;
			for (k=0; k<m->exon_start.size();k++){
				if (m->exon_start[k]<=m->genome_pos[j] && m->exon_end[k]>=m->genome_pos[j]){
					exon_id=k;
					break;
				}
			}
			if (exon_id!=-1){
				m->allele_exon[j]=exon_id;
				p->allele_exon[j]=exon_id;
				m->num_alleles_per_exon[exon_id]++;
				p->num_alleles_per_exon[exon_id]++;

			}
			else {
				fprintf(stdout,"error\n");
				exit(0);
			}
		}
	}
	return 0;
}

int Earrings::transcript_helper(string seq_filename,string gtf_file, 
							   vector<Transcript *> &transcripts){
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

int Earrings::align_reads(){
	int i,j,k;
	bool is_compatible=false;
	int num_maternal_alleles;
	int num_paternal_alleles;
	int total_alleles;
	set<BamAlignment *> unaligned;
	set<BamAlignment *> cleared;
	vector<set<BamAlignment *> >noninfo;
	vector<int> paternal_alleles;
	vector<int> maternal_alleles;

	noninfo.resize(maternal_transcripts.size());
	for(i=0;i<bam_reads.size();i++){
		is_compatible=false;
		for (j=0;j<maternal_transcripts.size();j++){
			// the maternal_transcript should be the same with the paternal_transcript
			if (maternal_transcripts[j]->is_compatible(bam_reads[i])){
				is_compatible=true;
				maternal_transcripts[j]->match_alleles(bam_reads[i],
													   total_alleles,
													   num_maternal_alleles,
													   maternal_alleles);
				paternal_transcripts[j]->match_alleles(bam_reads[i],
													   total_alleles,
													   num_paternal_alleles,
													   paternal_alleles);
				//fprintf(stdout,"%d\n",total_alleles);
				if ( total_alleles>0){
					if (num_paternal_alleles == num_maternal_alleles)
						noninfo[j].insert(bam_reads[i]);
					else {
						cleared.insert(bam_reads[i]);
						if (num_maternal_alleles!=0){
							maternal_transcripts[j]->register_read(bam_reads[i]);
						}
						else{
							paternal_transcripts[j]->register_read(bam_reads[i]);
						}
					}
				}
				else {
					noninfo[j].insert(bam_reads[i]);
				}
									
				//string mret=maternal_transcripts[j]->get_aligned_seq(bam_reads[i]);
				//string pret=paternal_transcripts[j]->get_aligned_seq(bam_reads[i]);
				//fprintf(stdout,"%s\n",ret.c_str());
			
	//fprintf(stdout,"%s\n",bam_reads[i]->QueryBases.c_str());
			}
		}
		if (is_compatible==false){
			unaligned.insert(bam_reads[i]);
			string c=get_cigar_string(*bam_reads[i]);
			if (c.find('N')!=c.npos &&bam_reads[i]->Name=="UNC9-SN296_0254:5:1201:10502:167862#TTAGGC"){
				maternal_transcripts[0]->output_segments();
				fprintf(stdout,"%s\n",(bam_reads[i]->Name.c_str()));
				fprintf(stdout,"%d\n",(bam_reads[i]->Position+1));
				fprintf(stdout,"%d\n",maternal_transcripts[0]->get_transcript_location(bam_reads[i]->Position+1));
				fprintf(stdout,"%s\n",get_cigar_string(*bam_reads[i]).c_str());
			}

			//fprintf(stdout,"%s\n",bam_reads[i]->QueryBases.c_str());
			//fprintf(stdout,"error\n");
		}
	}
	FILE *finfo=fopen(string(info->folder+"/"+info->gene_id+".landscape.plot.meta").c_str(),"w+");

	for (i=0;i<maternal_transcripts.size();i++){
		FILE *foutput;
		foutput=fopen(string(info->folder+"/"+maternal_transcripts[i]->name+".landscape.plot.info").c_str(),"w+");
		string name=maternal_transcripts[i]->name;
		LandscapePlot lp(maternal_transcripts[i],paternal_transcripts[i],noninfo[i]);
		lp.generate_landscape_plot(finfo, foutput);
		fclose(foutput);
	}
	fclose(finfo);
	//fprintf(stdout,"Unaligned\tUncleared\tCleared\tNoninfo\tTotal\n");
	//fprintf(stdout,"%d\t%d\t%d\t%d\t%d\n",unaligned.size(),
	//uncleared.size(),
	//cleared.size(),noninfo.size(),bam_reads.size());
	
	
}



int Earrings::test_allele_specific_transcript(){
	int i,j;
	for (i=0;i<maternal_transcripts.size();i++){
		bool maternal_dominance=false, paternal_dominance=false;
		vector<int>& maternal=maternal_transcripts[i]->num_info_reads_per_exon;
		vector<int>& paternal=paternal_transcripts[i]->num_info_reads_per_exon;
		for (j=0;j<maternal.size();j++){
			//if (maternal[j]>(paternal[j]+1)*3){
			if (maternal[j]>10 && paternal[j]==0){
				maternal_dominance=true;
			}
			//if (paternal[j]>(maternal[j]+1)*3){
			if (paternal[j]>10 && maternal[j]==0){
				paternal_dominance=true;
			} 
		}
		if (maternal_dominance || paternal_dominance) {
			fprintf(stdout,"here\n");
		}
		
	}

}


int Earrings::build_graph(){
	int i;
	Graph graph;
	FILE *foutput=fopen(string(info->folder+"/"+info->gene_id+".allele.specific.graph").c_str(),"w+");
	for (i=0;i<maternal_transcripts.size();i++){
		maternal_transcripts[i]->add_transcript_to_graph(&graph);
		paternal_transcripts[i]->add_transcript_to_graph(&graph);
	}
	graph.dump_graph(stdout);
	fclose(foutput);
}
