

#include "read_matcher.hpp"


void ReadMatcher::output(JewelerAlignment *al) {
	size_t i;
	if (mismatch_read_locations.size() == 0 ) return ;
	for ( i = 0; i < mismatch_read_locations.size(); i++) {
		fprintf(stdout, "%d\t", al->read_position[mismatch_read_locations[i]]);
	}
	for ( i = 0; i < mismatch_read_locations.size(); i++) {
		fprintf(stdout, "%d\t", mismatch_read_locations[i]);
	}
	for ( i = 0; i < mismatch_transcript_locations.size(); i++) {
		fprintf(stdout, "%d\t", mismatch_transcript_locations[i]);
	}
	fprintf(stdout, "\n");
}
