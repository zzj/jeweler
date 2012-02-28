#include "jeweler_alignment.hpp"

int get_read_position(JewelerAlignment *al, int i){
	if (al->IsReverseStrand()){
		return (al->QueryBases.size() - i -1);
	}
	else{
		return i;
	}
}
