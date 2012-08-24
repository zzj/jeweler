#include "transcript_mismatcher.hpp"

TranscriptMismatcher::TranscriptMismatcher() {

}

int TranscriptMismatcher::initialize() {
	size_t size = genome_pos2idx.size();
	int idx=0;
	for (auto i = genome_pos2idx.begin();
		 i != genome_pos2idx.end();
		 i++) {
		i->second=idx;
		idx++;
	}
	coverage.resize( size, 0);
	mismatches.resize( size, 0);
	mismatched_reads.resize(size);
	read_mismatch_qualities.resize(size);
	read_call_qualities.resize(size);
	read_mismatch_locations.resize(size);
	read_call_locations.resize(size);
	A_mismatches.resize( size, 0);
	C_mismatches.resize( size, 0);
	G_mismatches.resize( size, 0);
	T_mismatches.resize( size, 0);
	N_mismatches.resize( size, 0);
	return 0;
}

int TranscriptMismatcher::add_transcript(Transcript *t, int origin) {
	for (auto i=t->genome_pos.begin();
		 i != t->genome_pos.end(); 
		 i++) {
		genome_pos2idx[(*i)]=0;
		if (origin==TRANSCRIPT_PATERNAL) {
			genome_pos2paternal[(*i)]=t->seq[(i-t->genome_pos.begin())];
		}
		if (origin==TRANSCRIPT_MATERNAL) {
			genome_pos2maternal[(*i)]=t->seq[i-t->genome_pos.begin()];
		}
	}
	return 0;
}

int TranscriptMismatcher::add_mismatches(Transcript *transcript, JewelerAlignment *al, 
										 ReadMatcher *rm) {
	if (rm->mismatchars.size()> 10) {
		// TODO study why there are so many mismatches for the
		// alignment 
		return 0;
	}

	if (reads.find(al) == reads.end()) {
		for (size_t i = 0; i < rm->transcript_aligned_locations.size(); i++) {
			int idx = genome_pos2idx[transcript->genome_pos[rm->transcript_aligned_locations[i]]];

			coverage[ idx ]++;
			read_call_qualities[idx].push_back(al->Qualities[rm->read_aligned_locations[i]]);
			read_call_locations[idx].push_back(al->read_position[rm->read_aligned_locations[i]]);
		}
		reads.insert(al);
	}
	for (size_t i = 0; i < rm->mismatch_transcript_locations.size(); i++) {
		int idx = genome_pos2idx[transcript->genome_pos[rm->mismatch_transcript_locations[i]]];
		if (mismatched_reads[idx].find( al ) !=	mismatched_reads[idx].end()) {
			// have inserted before
			continue;
		}
		mismatched_reads[ idx ].insert( al );

		read_mismatch_qualities[ idx ].push_back( rm->mismatch_qualities[ i ] );
		read_mismatch_locations[ idx ].push_back( al->read_position[rm->mismatch_read_locations[ i ]] );

		char mismatchar = rm->mismatchars[i];
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
		case 'N': 
			N_mismatches[idx]++;
			break;
		default:
			fprintf(stderr, "%c is unknown sequence charactar\n", mismatchar);
		}
	}
	return 0;
}

int TranscriptMismatcher::dump(FILE *file) {
	size_t i;
	fprintf(file,"location\tcoverage\tmismatches\tA\tT\tC\tG\tN\tMaternal\tPaternal\n");
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++) {
		fprintf(file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\n",
				j->first, coverage[i], mismatches[i],
				A_mismatches[i],T_mismatches[i], C_mismatches[i], G_mismatches[i],
				N_mismatches[i],
				m->second, p->second);
		j++; m++; p++;
	}
	return 0;
}

int TranscriptMismatcher::write(FILE *file) {
	size_t i,k;
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++) {
		fprintf(file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c",
				j->first, coverage[i], mismatches[i],
				A_mismatches[i],T_mismatches[i], C_mismatches[i], G_mismatches[i],
				N_mismatches[i],
				m->second, p->second);
		for ( k = 0; k < read_mismatch_qualities[i].size(); k ++) {
			fprintf(file,"\t%c\t%d", read_mismatch_qualities[i][k],
					                 read_mismatch_locations[i][k]);
		}
		for ( k = 0; k < read_call_qualities[i].size(); k ++) {
			fprintf(file,"\t%c\t%d", read_call_qualities[i][k],
                                     read_call_locations[i][k]);
		}
		fprintf(file,"\n");
		j++; m++; p++;
	}
	return 0;
}

int TranscriptMismatcherAnalyzer::add_calls_by_quality(FILE * file, int num, 
													   map<char, int> & target_quality) {
	char tempc; int tempn;
	for (int i = 0; i < num; i ++ ) {
		fscanf(file, "\t%c\t%d", &tempc, &tempn);
		if (tempn<10) continue;
		if (target_quality.find(tempc) != target_quality.end()) {
			target_quality[tempc] ++;
		}
		else {
			target_quality[tempc] = 0 ;
		}
	}
	return 0;
}

void TranscriptMismatcherAnalyzer::append(FILE * file, string gene_id) {
	if (is_initialized) {
		fprintf(stdout, "TranscriptMismatcherAnalyzer: cannot load any new files\n");
		exit(0);
	}
	char maternal_char, paternal_char;
	int a,c,g,t,n;
	int coverage, genome_location, miss;
	map<char, int> tempv_quality;


	int ret = 0;
	while((ret = fscanf(file,"%d%d%d%d%d%d%d%d\t%c\t%c",
						&genome_location, &coverage, &miss,
						&a,&t, &c, &g, &n, 
						&maternal_char, &paternal_char) )== 10) {


		genome_locations.push_back( genome_location );
		coverages.push_back( coverage );
		mismatches.push_back( miss );
		if ( num_mismatches_histogram.find(miss) == num_mismatches_histogram.end() ) {
			num_mismatches_histogram[miss] = 0;
		}
		else {
			num_mismatches_histogram[miss] ++;
		}
		maternal_seq.push_back( maternal_char );
		paternal_seq.push_back( paternal_char );
		
		if (miss >= coverage * 0.02) {
			tempv_quality.clear();
			add_calls_by_quality(file, miss, tempv_quality);
			read_mismatch_qualities.push_back(tempv_quality);
			tempv_quality.clear();
			add_calls_by_quality(file, coverage, tempv_quality);
			read_qualities.push_back(tempv_quality);
		}
		else {
			skipped_locations.push_back(coverages.size()-1);
			add_calls_by_quality(file, miss, num_mismatches_by_quality_baseline);
			add_calls_by_quality(file, coverage, num_calls_by_quality_baseline);
			tempv_quality.clear();
			read_mismatch_qualities.push_back(tempv_quality);
			tempv_quality.clear();
			read_qualities.push_back(tempv_quality);
		}
		gene_ids.push_back(gene_id);

	}

}

TranscriptMismatcherAnalyzer::TranscriptMismatcherAnalyzer(string filename) {
	this -> filename = filename;
	
	is_initialized = false;
}
TranscriptMismatcherAnalyzer::TranscriptMismatcherAnalyzer(string filename, 
														   JewelerInfo *jeweler_info) {
	
	is_initialized = false;
	this -> filename = filename;
	for (size_t i=0; i < jeweler_info->gene_id.size(); i++) {
		string filename = string(jeweler_info->result_folder+"/"+jeweler_info->gene_id[i] +"/"+jeweler_info->gene_id[i]+".mismatcher.extra");
		FILE * finfo=fopen(filename.c_str(),"r");

		if (finfo == NULL) {
			fprintf (stdout, "cannot open file name %s\n", filename.c_str());
			continue;
		}
		append(finfo, jeweler_info->gene_id[i]);
		fclose(finfo);
		if ( i % 1 == 0) {
			fprintf(stdout, "%zu/%zu\n", i, jeweler_info->gene_id.size());
	 	}
	}
	end_loading();
}

void TranscriptMismatcherAnalyzer::end_loading() {
	is_initialized = true;
	p_values.resize(genome_locations.size(), 0);
	is_consistent_mismatches.resize(genome_locations.size(), false);
	is_skipped_locations.resize(genome_locations.size(), false);
	for (size_t i = 0; i <skipped_locations.size(); i ++) {
		is_skipped_locations.set(skipped_locations[i]);
	}
	// assuming paried reads are same
	// TODO: alow user to specify the length of the read
	// TODO: use input informatioin to get the number of reads
	num_reads = accumulate(coverages.begin(), coverages.end(), 0) / 100;
	num_locations = genome_locations.size();

}


// return the number of new found consistent misamtches
int TranscriptMismatcherAnalyzer::mark_consistent_mismatch() {
	int ret = 0;
	int i;
	for ( i = 0; i < num_locations; i ++) {
		if (is_skipped_locations.test(i)) {
			continue;
		}

		if (is_consistent_mismatches.test( i ) ) {
			continue;
		}
		
		if (p_values[i] < 0.00001) {
			is_consistent_mismatches.set( i );
			ret ++;
		}
	}
	return ret;
}

void TranscriptMismatcherAnalyzer::analyze() {
	int total = 0, ret;
	calculate_error<char>(read_qualities, read_mismatch_qualities, 
						  num_calls_by_quality, num_calls_by_quality_baseline, 
						  num_mismatches_by_quality, num_mismatches_by_quality_baseline,
						  error_rate_by_quality
						  );

	FILE * fd = fopen((filename+".quality.error.first").c_str(), "w+");
	dump_error_rate_by_quality (fd);
	fclose(fd);
	do {
		calculate_error<char>(read_qualities, read_mismatch_qualities, 
							  num_calls_by_quality, num_calls_by_quality_baseline, 
							  num_mismatches_by_quality, num_mismatches_by_quality_baseline,
							  error_rate_by_quality
							  );
		

		FILE * fd = fopen((filename+".quality.error.last").c_str(), "w+");
		dump_error_rate_by_quality(fd);
		fclose(fd);

		calculate_p_value<char>(read_qualities,error_rate_by_quality);
		ret = mark_consistent_mismatch() ;
		total += ret;
		fprintf(stdout, "there are %d of %d locations consistent mismatches\n", total, num_locations);
	} while ( ret >  0 );
	fd = fopen((filename+".locations").c_str(), "w+");
	dump_location_results (fd);
	fclose(fd);
	fd = fopen((filename+".consistent.locations").c_str(), "w+");
	dump_location_results (fd, true);
	fclose(fd);
	fd = fopen((filename+".log").c_str(), "w+");
	fprintf(fd, "%d\t%d\n", total, num_locations);
	fclose(fd);

}

void TranscriptMismatcherAnalyzer::dump_error_rate_by_quality(FILE * fd) {
	fprintf(fd, "phred\tnum.quality\tnum.mismatches\tnum.error\n");
	for ( auto i = num_calls_by_quality.begin(); 
		  i !=  num_calls_by_quality.end(); 
		  i ++) {
		fprintf(fd, "%c\t%d\t%d\t%e\n", 
				i->first, i->second, num_mismatches_by_quality[i->first], 
				error_rate_by_quality[ i->first]);
	}
}

void TranscriptMismatcherAnalyzer::dump_location_results(FILE *fd, bool only_yes) {
	int i;
	fprintf(fd, "location\tpvalue\tconsistent\n");
	for ( i = 0; i < num_locations; i ++) {
		if (only_yes && !is_consistent_mismatches.test(i)) continue;
		fprintf(fd, "%s\t%d\t%e\t%s\t%d\t%d\n", 
				gene_ids[i].c_str(),
				genome_locations[i], p_values[i], 
				is_consistent_mismatches.test(i) ? "Yes": "No", 
				coverages[i], mismatches[i]);

	}
}
