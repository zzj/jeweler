#include "transcript_mismatcher.hpp"

TranscriptMismatcher::TranscriptMismatcher(){
	
}


int TranscriptMismatcher::initialize(){
	size_t size = genome_pos2idx.size();
	int idx=0;
	for (auto i = genome_pos2idx.begin();
		 i != genome_pos2idx.end();
		 i++){
		i->second=idx;
		idx++;
	}
	coverage.resize(size,0);
	mismatches.resize( size, 0);
	mismatched_reads.resize(size);
	A_mismatches.resize( size, 0);
	C_mismatches.resize( size, 0);
	G_mismatches.resize( size, 0);
	T_mismatches.resize( size, 0);
	return 0;
}

int TranscriptMismatcher::add_transcript(Transcript *t, int origin){
	for (auto i=t->genome_pos.begin();
		 i != t->genome_pos.end(); 
		 i++){
		genome_pos2idx[(*i)]=0;
		if (origin==TRANSCRIPT_PATERNAL){
			genome_pos2paternal[(*i)]=t->seq[(i-t->genome_pos.begin())];
		}
		if (origin==TRANSCRIPT_MATERNAL){
			genome_pos2maternal[(*i)]=t->seq[i-t->genome_pos.begin()];
		}
	}
	return 0;
}

int TranscriptMismatcher::add_mismatches(Transcript *transcript, BamAlignment *al, 
										 vector<int> &locations,
										 vector<int> &mismatch_locations,
										 vector<char> & mismatchars){
	if (reads.find(al) == reads.end()){
		for (int i = 0; i < locations.size(); i++){
			coverage[genome_pos2idx[transcript->genome_pos[locations[i]]]]++;
		}
		reads.insert(al);
	}
	for (int i = 0; i < mismatch_locations.size(); i++){
		int idx = genome_pos2idx[transcript->genome_pos[mismatch_locations[i]]];

		if (mismatched_reads[idx].find( al ) !=	mismatched_reads[idx].end()){
			// have inserted before
			continue;
		}
		mismatched_reads[idx].insert(al);

		char mismatchar = mismatchars[i];
		mismatches[idx]++;

		switch (mismatchar) {
		case 'A': 
			A_mismatches[idx]++;
			break;
		case 'T': 
			T_mismatches[idx]++;
			break;
		case 'C': 
			C_mismatches[idx]++;
			break;
		case 'G': 
			G_mismatches[idx]++;
			break;
		default:
			fprintf(stderr, "%c is unknown sequence charactar\n", mismatchar);
		}
	}
	return 0;
}


int TranscriptMismatcher::dump(FILE *file){
	int i;
	fprintf(file,"locations\ttotal\tmismatches\tA\tT\tC\tG\n");
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++){
		fprintf(file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\n",
				j->first, coverage[i], mismatches[i],
				A_mismatches[i],T_mismatches[i], C_mismatches[i], G_mismatches[i],
				m->second, p->second);
		j++; m++; p++;
	}
	return 0;
}
