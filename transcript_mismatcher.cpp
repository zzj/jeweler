#include "transcript_mismatcher.hpp"
#include "jeweler_info.hpp"
#include "transcript.hpp"
#include "constants.hpp"
#include "common.hpp"
#include "zleveldb.hpp"

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
	coverage.resize(size, 0);
	mismatches.resize(size, 0);
	mismatched_reads.resize(size);
	read_mismatch_qualities.resize(size);
	read_call_qualities.resize(size);
	read_mismatch_locations.resize(size);
	read_call_locations.resize(size);
	A_mismatches.resize(size, 0);
	C_mismatches.resize(size, 0);
	G_mismatches.resize(size, 0);
	T_mismatches.resize(size, 0);
	N_mismatches.resize(size, 0);
	return 0;
}

int TranscriptMismatcher::add_transcript(Transcript *t, int origin) {
	for (auto i=t->genome_pos.begin();
		 i != t->genome_pos.end();
		 i++) {
		genome_pos2idx[(*i)]=0;
		if (origin == TRANSCRIPT_PATERNAL) {
			genome_pos2paternal[(*i)] = t->seq()[(i-t->genome_pos.begin())];
		}
		if (origin == TRANSCRIPT_MATERNAL) {
			genome_pos2maternal[(*i)] = t->seq()[i-t->genome_pos.begin()];
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
			read_call_locations[idx].push_back(rm->read_aligned_locations[i]);
		}
        reads.insert(al);
    }
    for (size_t i = 0; i < rm->mismatch_transcript_locations.size(); i++) {
        int idx = genome_pos2idx[transcript->genome_pos[rm->mismatch_transcript_locations[i]]];
        if (mismatched_reads[idx].find(al) != mismatched_reads[idx].end()) {
            // have inserted before
            continue;
        }
        mismatched_reads[ idx ].insert(al);

        read_mismatch_qualities[ idx ].push_back(rm->mismatch_qualities[i]);
        read_mismatch_locations[ idx ].push_back(rm->mismatch_read_locations[i]);

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

void TranscriptMismatcher::dump(Jeweler::EarringsData::Mismatcher *mismatcher) {
	size_t i,k;
    string paternal_chars;
    string maternal_chars;
	auto j=genome_pos2idx.begin();
	auto m=genome_pos2maternal.begin();
	auto p=genome_pos2paternal.begin();
	for (i = 0 ; i < mismatches.size() ; i++) {
        Jeweler::EarringsData::Mismatcher::Mismatch *data = 
            mismatcher->add_mismatch();
        data->set_genome_position(j->first);
        data->set_coverage(coverage[i]);
        data->set_num_mismatches(mismatches[i]);
        data->set_num_a(A_mismatches[i]);
        data->set_num_t(T_mismatches[i]);
        data->set_num_c(C_mismatches[i]);
        data->set_num_g(G_mismatches[i]);
        data->set_paternal_char(string(1,m->second));
        data->set_maternal_char(string(1,m->second));
        string quality_string;
        for (k = 0; k < read_mismatch_qualities[i].size(); k ++) {

            Jeweler::EarringsData::Mismatcher::Mismatch::Call *call = data->add_call();
            quality_string += read_mismatch_qualities[i][k];
            call->set_read_position(read_mismatch_locations[i][k]);
            call->set_is_mismatch(true);
		}
        for (k = 0; k < read_call_qualities[i].size(); k ++) {
            Jeweler::EarringsData::Mismatcher::Mismatch::Call *call = data->add_call();
            quality_string += read_call_qualities[i][k];
            call->set_read_position(read_call_locations[i][k]);
            call->set_is_mismatch(false);
		}
        data->set_quality_string(quality_string);
		j++; m++; p++;
	}



	return ;
}

int TranscriptMismatcherAnalyzer::add_calls_by_quality(
                   const Jeweler::EarringsData::Mismatcher::Mismatch &data,
                   map<char, int> & target_mismatch_quality,
                   map<char, int> & target_match_quality) {

    assert(data.call_size() == data.quality_string().size());

	for (int i = 0; i < data.call_size(); i ++) {
        // ignore the first 10 characters.
		if (data.call(i).read_position() < 10) continue;
        if (data.call(i).is_mismatch()) {
            map_add_count<char>(target_mismatch_quality, data.quality_string()[i]);
        }
        else {
            map_add_count<char>(target_match_quality, data.quality_string()[i]);
        }
	}
	return 0;
}

void TranscriptMismatcherAnalyzer::append(const Jeweler::EarringsData::Mismatcher &data,
                                          string gene_id) {
	if (is_initialized) {
		fprintf(stdout, "TranscriptMismatcherAnalyzer: cannot load any new files\n");
		exit(0);
	}

    for (int i = 0; i < data.mismatch_size(); i++) {
        map<char, int> tempv_mis_quality;
        map<char, int> tempv_match_quality;
        const Jeweler::EarringsData::Mismatcher::Mismatch &mismatch = 
            data.mismatch(i);
        int miss = mismatch.num_mismatches();
        int coverage = mismatch.coverage();
		genome_locations.push_back(mismatch.genome_position());
		coverages.push_back(coverage);
		mismatches.push_back(miss);
		if (num_mismatches_histogram.find(miss) == num_mismatches_histogram.end()) {
			num_mismatches_histogram[miss] = 0;
		}
		else {
			num_mismatches_histogram[miss] ++;
		}
		maternal_seq.push_back(mismatch.maternal_char()[0]);
		paternal_seq.push_back(mismatch.paternal_char()[0]);

        add_calls_by_quality(mismatch, tempv_mis_quality, tempv_match_quality);

		if (miss >= coverage * 0.02) {
			read_mismatch_qualities.push_back(tempv_mis_quality);
			read_qualities.push_back(tempv_match_quality);
		}
		else {
			skipped_locations.push_back(coverages.size()-1);
			add_calls_by_quality(mismatch,
                                 num_mismatches_by_quality_baseline,
                                 num_calls_by_quality_baseline);
			tempv_mis_quality.clear();
			read_mismatch_qualities.push_back(tempv_mis_quality);
			tempv_match_quality.clear();
			read_qualities.push_back(tempv_match_quality);
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
    shared_ptr<ZMegaFile> zmf (new ZMegaFile(jeweler_info->result_folder));
	is_initialized = false;
	this -> filename = filename;
	for (size_t i=0; i < jeweler_info->gene_id.size(); i++) {
        shared_ptr<Jeweler::EarringsData> ed = 
            zmf->get<Jeweler::EarringsData>(jeweler_info->gene_id[i]);
        if (ed.get() == NULL) {
            continue;
        }
		append(ed->mismatcher(), jeweler_info->gene_id[i]);
		if (i % 100 == 0) {
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
	for (i = 0; i < num_locations; i ++) {
		if (is_skipped_locations.test(i)) {
			continue;
		}

		if (is_consistent_mismatches.test(i)) {
			continue;
		}

		if (p_values[i] < 0.00001) {
			is_consistent_mismatches.set(i);
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
	} while (ret >  0);
	fd = fopen((filename+".locations").c_str(), "w+");
	dump_location_results (fd);
	fclose(fd);
	fd = fopen((filename+".consistent.locations").c_str(), "w+");
	dump_location_results (fd, true);
	fclose(fd);
}

void TranscriptMismatcherAnalyzer::dump_error_rate_by_quality(FILE * fd) {
	fprintf(fd, "phred\tnum.quality\tnum.mismatches\tnum.error\n");
	for (auto i = num_calls_by_quality.begin();
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
	for (i = 0; i < num_locations; i ++) {
		if (only_yes && !is_consistent_mismatches.test(i)) continue;
		fprintf(fd, "%s\t%d\t%e\t%s\t%d\t%d\n",
				gene_ids[i].c_str(),
				genome_locations[i], p_values[i],
				is_consistent_mismatches.test(i) ? "Yes": "No",
				coverages[i], mismatches[i]);

	}
}
