#ifndef _PROBABILITY_HPP_
#define _PROBABILITY_HPP_

#include <vector> 
#include <cmath> 
#include <cstdio>
using namespace std;



// calcuate the prod of probabilites of p
// log(\prod p(x_i))  =  \sum( log p(x_i))

double log_prod_prob(vector<double> &prod, bool is_complementary = false);

// calcuate the prod of probabilites of p
// log(\sum_{i} \prod_{j} (p(x_j)I(i != j) + (1-p(x_j))I(i == j)))
double log_prod_prob_with_one_flipped(vector<double> &prod, bool is_complementary = false);


#endif
