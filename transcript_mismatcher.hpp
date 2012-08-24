#ifndef _TRANSCRIPT_MISMATCHER_H_
#define _TRANSCRIPT_MISMATCHER_H_

#include "transcript.hpp"
#include "constants.hpp"
#include "jeweler.hpp"
#include <boost/dynamic_bitset.hpp>
#include <boost/math/distributions/poisson.hpp>

#include "math/probability.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>

using namespace std;
class JewelerInfo;

class TranscriptMismatcher{

public:

	TranscriptMismatcher();

	int initialize();

	int add_transcript(Transcript *, int origin);

	int add_mismatches(Transcript *transcript, JewelerAlignment *al,
					   ReadMatcher *rm);

	int dump(FILE *);
	int write(FILE *);

	map<int, int> genome_pos2idx;
	map<int, char> genome_pos2paternal;
	map<int, char> genome_pos2maternal;

	vector<int> mismatches; // number of non informative mismatches
							// per base
	set<JewelerAlignment *> reads;
	vector<int> coverage;
	// the list of mismatching reads per base pair
	vector<set<JewelerAlignment *> > mismatched_reads;
	// the list of mismatching qualities
	vector<vector<char > > read_mismatch_qualities;
	vector<vector<char > > read_call_qualities;

	// the list of mismatching locations
	vector<vector<int > > read_mismatch_locations;
	vector<vector<int > > read_call_locations;

	// number of non informative mismatches equaling to A per base
	vector<int> A_mismatches;
	// number of non informative mismatches equaling to C per base
	vector<int> C_mismatches;
	// number of non informative mismatches equaling to G per base
	vector<int> G_mismatches;
	// number of non informative mismatches equaling to T per base
	vector<int> T_mismatches;
	// number of non informative mismatches equaling to N per base
	vector<int> N_mismatches;

};

class TranscriptMismatcherAnalyzer{
public:
	vector<int> genome_locations;
	vector<int> coverages;
	vector<int> mismatches;
	vector<char> maternal_seq, paternal_seq;
	vector<map<char, int> > read_mismatch_qualities;
	vector<map<char, int> > read_qualities;
	vector<string> gene_ids;
	map<char, double> error_rate_by_quality;
	// the list of mismatching locations
	//vector<vector<int > > read_mismatch_locations;
	//vector<vector<int > > read_locations;
	//// map<int, double> error_rate_by_location;
	// // // number of non informative mismatches equaling to A per base
	// // vector<int> A_mismatches;
	// // // number of non informative mismatches equaling to C per base
	// // vector<int> C_mismatches;
	// // // number of non informative mismatches equaling to G per base
	// // vector<int> G_mismatches;
	// // // number of non informative mismatches equaling to T per base
	// // vector<int> T_mismatches;
	// // // number of non informative mismatches equaling to N per base
	// // vector<int> N_mismatches;


	// // prefix of output files
	string filename;
	int num_reads;
	int num_locations;
	boost::dynamic_bitset<> is_consistent_mismatches;

	// if the location has less than 1% mismatches. put the count into
	// the baseline, and skip the current location.
	boost::dynamic_bitset<> is_skipped_locations;
	vector<int> skipped_locations;
	map<char, int> num_mismatches_by_quality;
	map<char, int> num_calls_by_quality;

	map<char, int> num_mismatches_by_quality_baseline;
	map<char, int> num_calls_by_quality_baseline;

	map<int, int> num_mismatches_by_location;
	map<int, int> num_calls_by_location;
	map<int, int> num_mismatches_histogram;
	vector<double> p_values;

	bool is_initialized;

	TranscriptMismatcherAnalyzer(string filename);
	TranscriptMismatcherAnalyzer(string filename, JewelerInfo *jeweler_info);

	void append(FILE *, string gene_id );

	void end_loading();

    template<class T>
    void calculate_error(const vector<map<T, int> > &read_calls,
                        const vector<map<T, int> > &read_mismatches,
                        map<T, int> &num_calls,
                        map<T, int> &num_calls_baseline,
                        map<T, int> &num_mismatches,
                        map<T, int> &num_mismatches_baseline,
                        map<T,double> &error_rate
                        );
    int add_calls_by_quality(FILE * file, int num,
                             map<char, int> & target_quality);

    template<class T>
    void calculate_p_value(const vector<map<T, int> > &read_calls,
                          map<T,double> &error_rate
                          );

    int mark_consistent_mismatch();

    void dump_error_rate_by_quality(FILE *);
    void dump_location_results(FILE *, bool only_yes = false);

    void analyze();
};

template<class T>
void TranscriptMismatcherAnalyzer::calculate_p_value(const vector<map<T, int> > &read_calls,
													map<T,double> &error_rate
													){
	int i;
	double le_cam_upper_bound;
	double mean;
	vector<double> prob;

	for ( i = 0; i < num_locations; i ++) {
		prob.clear();
		if ( is_consistent_mismatches.test( i ) ){
			continue;
		}
		le_cam_upper_bound = 0;
		mean = 0;
		for (auto j = read_calls[i].begin() ;
			 j != read_calls[i].end();
			 j ++){
			le_cam_upper_bound = pow(error_rate[ j->first ],2);
			mean += error_rate[j->first ];
			for ( int k = 0; k < j->second; k++){
				prob.push_back(error_rate[ j->first ]);
			}
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
template<class T>
void TranscriptMismatcherAnalyzer::calculate_error(const vector<map<T, int> > &read_calls,
												  const vector<map<T, int> > &read_mismatches,
												  map<T, int> &num_calls,
												  map<T, int> &num_calls_baseline,
												  map<T, int> &num_mismatches,
												  map<T, int> &num_mismatches_baseline,
												  map<T,double> &error_rate
												  ){
	int i;
	// num of mismatches by quality in read
	num_mismatches = num_mismatches_baseline;
	num_calls = num_calls_baseline;
	error_rate.clear();

	for ( i = 0; i < num_locations; i ++){
		// Right now, just remove locations
		// TODO: remove the whole transcript if it contains consistent
		// mismatches
		if ( is_consistent_mismatches.test( i ) ){
			continue;
		}
		if (is_skipped_locations.test(i)){
			continue;
		}
		for (auto j = read_mismatches[i].begin() ;
			 j != read_mismatches[i].end();
			 j ++){
			if (num_mismatches.find(j->first) ==
				num_mismatches.end()){
				num_mismatches[ j->first ] = j->second;
			}
			else {
				num_mismatches[ j->first ] += j->second;
			}
		}
		for (auto j = read_calls[i].begin() ;
			 j != read_calls[i].end();
			 j ++){
			if (num_calls.find(j->first) ==
				num_calls.end()){
				num_calls[ j->first ] = j->second;
			}
			else {
				num_calls[j->first ] += j->second;
			}
		}
 	}

	for ( auto i = num_calls.begin(); i !=  num_calls.end(); i ++){
		if (num_mismatches.find(i->first) == num_mismatches.end()) {
			num_mismatches[ i->first ] = 0;
		}
		error_rate[ i->first ] = (double)  num_mismatches[i->first] / i->second ;
	}
}

#endif
