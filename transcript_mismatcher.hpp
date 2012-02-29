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
	vector<vector<char> > read_mismatch_qualities;
	vector<vector<char> > read_qualities;
	map<char, double> error_rate_by_quality;	
	// the list of mismatching locations 
	vector<vector<int > > read_mismatch_locations;
	vector<vector<int > > read_locations;
	map<int, double> error_rate_by_location;	
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
	

	// prefix of output files
	string filename;
	int num_reads;
	int num_locations;
	boost::dynamic_bitset<> is_consistent_mismatches;

	map<char, int> num_mismatches_by_quality;
	map<char, int> num_calls_by_quality;
	map<int, int> num_mismatches_by_location;
	map<int, int> num_calls_by_location;
	map<int, int> num_mismatches_histogram;
	vector<double> p_values;

	bool is_initialized;
	
	TranscriptMismatcherAnalyzer(string filename);
	TranscriptMismatcherAnalyzer(string filename, JewelerInfo *jeweler_info);

	int append(FILE *);
	
	int end_loading();
	
	template<class T>
	int calculate_error(const vector<vector<T> > &read_calls,
						const vector<vector<T> > &read_mismatches,
						map<T, int> &num_calls,
						map<T, int> &num_mismatches,
						map<T,double> &error_rate
						);

	template<class T> 
	int calculate_p_value(const vector<vector<T> > &read_calls,
						  map<T,double> &error_rate
						  );


	int mark_consistent_mismatch();

	int dump_error_rate_by_quality(FILE *);
	int dump_error_rate_by_location(FILE *);

	int dump_location_results(FILE *, bool only_yes = false);

	int analyze();
};

template<class T> 
int TranscriptMismatcherAnalyzer::calculate_p_value(const vector<vector<T> > &read_calls,
													map<T,double> &error_rate
													){
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
		for (j = 0 ; j < read_calls[i].size(); j ++){
			le_cam_upper_bound = pow(error_rate[ read_calls[ i ][ j ] ],2);
			mean += error_rate[ read_calls[ i ][ j ] ];
			prob.push_back(error_rate[ read_calls[ i ][ j ] ]);
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
int TranscriptMismatcherAnalyzer::calculate_error(const vector<vector<T> > &read_calls,
												  const vector<vector<T> > &read_mismatches,
												  map<T, int> &num_calls,
												  map<T, int> &num_mismatches,
												  map<T,double> &error_rate
												  ){
	int i, j;
	// num of mismatches by quality in read
	num_mismatches.clear();
	num_calls.clear();
	error_rate.clear();

	for ( i = 0; i < num_locations; i ++){
		if ( is_consistent_mismatches.test( i ) ){
			continue;
		}
		// Right now, just remove locations 
		// TODO: remove the whole transcript if it contains consistent
		// mismatches
		for (j = 0 ; j < read_mismatches[i].size(); j ++){
			if (num_mismatches.find(read_mismatches[ i ][ j ]) == 
				num_mismatches.end()){
				num_mismatches[ read_mismatches[ i ][ j ] ] = 1;
			}
			else {
				num_mismatches[ read_mismatches[ i ][ j ] ] ++;
			}
		}
		for (j = 0 ; j < read_calls[i].size(); j ++){
			if (num_calls.find(read_calls[ i ][ j ]) == 
				num_calls.end()){
				num_calls[ read_calls[ i ][ j ] ] = 1;
			}
			else {
				num_calls[ read_calls[ i ][ j ] ] ++;
			}
		}
 	}

	for ( auto i = num_calls.begin(); i !=  num_calls.end(); i ++){
		if (num_mismatches.find(i->first) == num_mismatches.end()) {
			num_mismatches[ i->first ] = 0;
		}
		error_rate[ i->first ] = (double)  num_mismatches[i->first] / i->second ;
	}
	return 0;
}

#endif
