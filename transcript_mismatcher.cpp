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
	coverage.resize( size, 0);
	mismatches.resize( size, 0);
	mismatched_reads.resize(size);
	read_mismatch_qualities.resize(size);
	read_qualities.resize(size);
	A_mismatches.resize( size, 0);
	C_mismatches.resize( size, 0);
	G_mismatches.resize( size, 0);
	T_mismatches.resize( size, 0);
	N_mismatches.resize( size, 0);
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

int TranscriptMismatcher::add_mismatches(Transcript *transcript, JewelerAlignment *al, 
										 vector<int> &transcript_locations,
										 vector<int> &read_locations,
										 vector<int> &mismatch_locations,
										 vector<char> &read_mismatch_qualities_in,
										 vector<char> & mismatchars){
	if (mismatchars.size()> 4) {
		// TODO study why there are so many mismatches for the
		// alignment 
		return 0;
	}

	if (reads.find(al) == reads.end()){
		for (int i = 0; i < transcript_locations.size(); i++){
			int idx = genome_pos2idx[transcript->genome_pos[transcript_locations[i]]];
			coverage[ idx ]++;
			read_qualities[ idx ].push_back( al->Qualities[ read_locations[i] ] );
		}
		reads.insert(al);
	}
	for (int i = 0; i < mismatch_locations.size(); i++){
		

		int idx = genome_pos2idx[transcript->genome_pos[mismatch_locations[i]]];

		if (mismatched_reads[idx].find( al ) !=	mismatched_reads[idx].end()){
			// have inserted before
			continue;
		}
		mismatched_reads[ idx ].insert( al );
		read_mismatch_qualities[ idx ].push_back( read_mismatch_qualities_in[ i ] );
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
		case 'N': 
			N_mismatches[idx]++;
			break;
		default:
			fprintf(stderr, "%c is unknown sequence charactar\n", mismatchar);
		}
	}
	return 0;
}

int TranscriptMismatcher::dump(FILE *file){
	int i;
	fprintf(file,"location\tcoverage\tmismatches\tA\tT\tC\tG\tN\tMaternal\tPaternal\n");
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++){
		fprintf(file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\n",
				j->first, coverage[i], mismatches[i],
				A_mismatches[i],T_mismatches[i], C_mismatches[i], G_mismatches[i],
				N_mismatches[i],
				m->second, p->second);
		j++; m++; p++;
	}
	return 0;
}

int TranscriptMismatcher::write(FILE *file){
	int i,k;
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++){
		fprintf(file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c",
				j->first, coverage[i], mismatches[i],
				A_mismatches[i],T_mismatches[i], C_mismatches[i], G_mismatches[i],
				N_mismatches[i],
				m->second, p->second);
		for ( k = 0; k < read_mismatch_qualities[i].size(); k ++){
			fprintf(file,"\t%c", read_mismatch_qualities[i][k]);
		}
		for ( k = 0; k < read_qualities[i].size(); k ++){
			fprintf(file,"\t%c", read_qualities[i][k]);
		}
		fprintf(file,"\n");
		j++; m++; p++;
	}
	return 0;
}

int TranscriptMismatcherAnalyzer::append(FILE * file){
	if (is_initialized){
		fprintf(stdout, "TranscriptMismatcherAnalyzer: cannot load any new files\n");
		exit(0);
	}
	char maternal_char, paternal_char;
	int a,c,g,t,n;
	int coverage, genome_location, miss;
	vector<char> tempv;
	int temp;
	int max=error_rate.size();
	int ret = 0;
	while((ret = fscanf(file,"%d%d%d%d%d%d%d%d\t%c\t%c",
						&genome_location, &coverage, &miss,
						&a,&t, &c, &g, &n, 
						&maternal_char, &paternal_char) )== 10){

		tempv.clear();
		genome_locations.push_back( genome_location );
		coverages.push_back( coverage );
		mismatches.push_back( miss );
		if ( num_mismatches_histogram.find(miss) == num_mismatches_histogram.end() ){
			num_mismatches_histogram[miss] = 0;
		}
		else {
			num_mismatches_histogram[miss] ++;
		}
		A_mismatches.push_back( a );
		T_mismatches.push_back( t );
		C_mismatches.push_back( c );
		G_mismatches.push_back( g );
		N_mismatches.push_back( n );
		maternal_seq.push_back( maternal_char );
		paternal_seq.push_back( paternal_char );
		
		for (int i = 0; i < miss; i ++ ){
			fscanf(file, "\t%c", &temp);
			tempv.push_back(temp);
		}
		read_mismatch_qualities.push_back(tempv);
		tempv.clear();
		
		for (int i = 0; i < coverage; i ++ ){
			fscanf(file, "\t%c", &temp);
			tempv.push_back(temp);
		}
		
		read_qualities.push_back(tempv);
	}

}

TranscriptMismatcherAnalyzer::TranscriptMismatcherAnalyzer(string filename){
	this -> filename = filename;
	is_initialized = false;
}
TranscriptMismatcherAnalyzer::TranscriptMismatcherAnalyzer(string filename, 
														   vector<TranscriptInfo *> &ti){
	
	is_initialized = false;
	this -> filename = filename;
	for ( int i=0; i < ti.size(); i++){
		TranscriptInfo * info = ti[ i ];
		string filename = string(info->folder+"/"+info->gene_id+".mismatcher.extra");
		FILE * finfo=fopen(filename.c_str(),"r+");

		if (finfo == NULL){
			fprintf (stdout, "cannot open file name %s\n", filename.c_str());
			exit(0);
		}
		append(finfo);
		fclose(finfo);
		if ( i % 100 == 0) {
			fprintf(stdout, "%d/%d\n", i, ti.size());
		}
	}
	end_loading();
}

int TranscriptMismatcherAnalyzer::end_loading(){
	is_initialized = true;
	p_values.resize(genome_locations.size(), 0);
	is_consistent_mismatches.resize(genome_locations.size(), false);
	// assuming paried reads are same
	// TODO: alow user to specify the length of the read
	// TODO: use input informatioin to get the number of reads
	num_reads = accumulate(coverages.begin(), coverages.end(), 0) / 100;
	num_locations = genome_locations.size();

}

int TranscriptMismatcherAnalyzer::calculate_error(){
	int i, j;
	// num of mismatches by location in read
	num_mismatches.clear();
	num_qualities.clear();
	error_rate.clear();

	for ( i = 0; i < num_locations; i ++){
		if ( is_consistent_mismatches.test( i ) ){
			continue;
		}
		// Right now, just remove locations 
		// TODO: remove the whole transcript if it contains consistent
		// mismatches
		for (j = 0 ; j < read_mismatch_qualities[i].size(); j ++){
			if (num_mismatches.find(read_mismatch_qualities[ i ][ j ]) == 
				num_mismatches.end()){
				num_mismatches[ read_mismatch_qualities[ i ][ j ] ] = 1;
			}
			else {
				num_mismatches[ read_mismatch_qualities[ i ][ j ] ] ++;
			}
		}
		for (j = 0 ; j < read_qualities[i].size(); j ++){
			if (num_qualities.find(read_qualities[ i ][ j ]) == 
				num_qualities.end()){
				num_qualities[ read_qualities[ i ][ j ] ] = 1;
			}
			else {
				num_qualities[ read_qualities[ i ][ j ] ] ++;
			}
		}
 	}

	for ( auto i = num_qualities.begin(); 
		  i !=  num_qualities.end(); 
		  i ++){
		if (num_mismatches.find(i->first) == num_mismatches.end()) {
			num_mismatches[ i->first ] = 0;
		}
		error_rate[ i->first ] = (double)  num_mismatches[i->first] / i->second ;
		fprintf(stdout, "%c\t%d\t%d\t%e\n", 
				i->first, i->second, num_mismatches[i->first], 
				error_rate[ i->first]);

	}
	return 0;
}

int TranscriptMismatcherAnalyzer::calculate_p_value(){
	int i;
	double le_cam_upper_bound;
	double mean;
	vector<double> prob;
	int j;

	for ( i = 0; i < num_locations; i ++) {
		prob.clear();
		if ( is_consistent_mismatches.test( i ) ){
			continue;
		}
		le_cam_upper_bound = 0;
		mean = 0;
		for (j = 0 ; j < read_qualities[i].size(); j ++){
			le_cam_upper_bound = pow(error_rate[ read_qualities[ i ][ j ] ],2);
			mean += error_rate[ read_qualities[ i ][ j ] ];
			prob.push_back(error_rate[ read_qualities[ i ][ j ] ]);
		}
		p_values[i] = 0;
		if (coverages[i] == 0){
			p_values[i] = 1;
		}

		if (mean == 0) continue;
		boost::math::poisson_distribution<double> pd = boost::math::poisson_distribution<double> (mean);
		bool test = false;
		if (coverages[i] == 0){
			p_values[i] = 1;
		}
		else if (mismatches[i] == 0){
			p_values[i] = 1;
		}
		else if (mismatches[i] == 1){
			p_values[i] = 1 - exp(log_prod_prob(prob, true));
		}
		else if (mismatches[i] ==2){
			p_values[i] = 1 - exp(log_prod_prob(prob, true)) - exp(log_prod_prob_with_one_flipped(prob, true));
		}
		else if (mismatches[i] == coverages[i]){
			p_values[i] = exp(log_prod_prob(prob));
		}
		else if (mismatches[i] == coverages[i] - 1){
			p_values[i] = exp(log_prod_prob(prob) )
				+ exp(log_prod_prob_with_one_flipped(prob));
		}
		else {
			test = true;
			p_values[i]=(cdf(complement(pd, mismatches[i] - 1)));
		}
		if (false){
			if (mismatches[i] == 0) continue;
			double test1 = p_values[i];
			double test2 = (cdf(complement(pd, mismatches[i] - 1)));
			if (abs(test1 -test2) > 0.00001){
				fprintf(stdout, "Cov %d Miss %d Exact %lf estimate %lf\n",
						coverages[i], mismatches[i], test1, test2);
			}
		}
		//if (le_cam_upper_bound > 0) fprintf(stdout, "%lf\n", le_cam_upper_bound);
	}
}

// return the number of new found consistent misamtches
int TranscriptMismatcherAnalyzer::mark_consistent_mismatch(){
	int ret = 0;
	int i;
	for ( i = 0; i < num_locations; i ++){
		if (is_consistent_mismatches.test( i ) ){
			continue;
		}
		if (p_values[i] < 0.00001){
			is_consistent_mismatches.set( i );
			ret ++;
		}
	}
	return ret;
}

int TranscriptMismatcherAnalyzer::analyze(){
	int total = 0, ret;
	calculate_error();

	FILE * fd = fopen((filename+".error.first").c_str(), "w+");
	dump_error_rate (fd);
	fclose(fd);
	do {
		calculate_error();
		FILE * fd = fopen((filename+".error.last").c_str(), "w+");
		dump_error_rate (fd);
		fclose(fd);
		calculate_p_value();
		ret = mark_consistent_mismatch() ;
		total += ret;
		fprintf(stdout, "there are %d of %d locations consistent mismatches\n", total, num_locations);
	} while ( ret >  0 );
	fd = fopen((filename+".locations").c_str(), "w+");
	dump_location_results (fd);
	fclose(fd);
	fd = fopen((filename+".log").c_str(), "w+");
	fprintf(fd, "%d\t%d\n", total, num_locations);
	fclose(fd);

}

int TranscriptMismatcherAnalyzer::dump_error_rate(FILE * fd){
	fprintf(fd, "phred\tnum.quality\tnum.mismatches\tnum.error\n");
	for ( auto i = num_qualities.begin(); 
		  i !=  num_qualities.end(); 
		  i ++){
		fprintf(fd, "%c\t%d\t%d\t%e\n", 
				i->first, i->second, num_mismatches[i->first], 
				error_rate[ i->first]);
	}
	return 0;
}

int TranscriptMismatcherAnalyzer::dump_location_results(FILE *fd){
	int i;
	fprintf(fd, "location\tpvalue\tconsistent\n");
	for ( i = 0; i < num_locations; i ++){
		fprintf(fd, "%d\t%e\t%s\t%d\t%d\n", genome_locations[i], p_values[i], 
				is_consistent_mismatches.test(i) ? "Yes": "No", coverages[i], mismatches[i]);
	}
}
