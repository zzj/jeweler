#include "probability.hpp"

// log(1-x)

double log1mx(double x){
	if (x < 0 || x > 1){
		fprintf(stdout, " out the range " );
	}
	if (x>0.00001){
		return log(1-x);
	}
	else {
		int i;
		double ret = 0;
		// use taylor series
		// log(1-x) = -\sum_{n}\frac{x^n}{n}

		for (i = 1; i < 20; i ++){
			ret += exp(x * i) / i;
		}
	}
}

double log_prod_prob(vector<double> &prod, bool is_complementary ){
	int n = prod.size();
	int i;
	double ret = 0;
	for ( i = 0; i< prod.size(); i ++){
		if (is_complementary) {
			ret += log1mx( prod[i] );
		}
		else {
			ret += log( prod[i] );
		}
	}
	return ret;
}

double log_prod_prob_with_one_flipped(vector<double> &prod, bool is_complememtary){

	// left[i] \sum_{j}^{i} log (prod[j])
	int i;
	double left = 0;
	double right = 0;
	double ret = 0;
	for ( i = 0; i < prod.size(); i ++){
		if (is_complememtary)
			right += log1mx(prod[i]);
		else {
			right += log(prod[i]);
		}
	}
	
	for ( i = 0; i < prod.size(); i ++){
		
		double curr = 0;
		if ( is_complememtary) {
			right -= log1mx(prod[i]);
			curr = left + right + log(prod[i]);
			left += log1mx(prod[i]);
		}
		else {
			right -= log(prod[i]);
			curr = left + right + log1mx(prod[i]);
			left += log(prod[i]);
		}
		
		ret += exp(curr);		
	}
	return log(ret);
}
