#include <RcppEigen.h>
#include "multinomial.h"
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
IntegerVector multinomial_resampling_n_(const NumericVector & weights, int ndraws){
  RNGScope scope;
  int nparticles = weights.size();
  // create vector of ancestor variables
  IntegerVector ancestors(ndraws);
  // cumulative weights
  NumericVector cumsumw = cumsum(weights);
  NumericVector uniforms = runif(ndraws);
  // the following avoids sorting the uniforms
  // it is based on the directed generation of a vector of sorted uniforms 
  // using scaled Exponential variables
  double sumw = cumsumw(nparticles - 1);
  double lnMax = 0;
  int j = nparticles;
  for (int i = ndraws; i > 0; i--){
    lnMax += log(uniforms(i-1)) / i;
    uniforms(i-1) = sumw * exp(lnMax);
    while (j > 0 && uniforms(i-1) < cumsumw(j-1)){
      j --;
    }
    ancestors(i-1) = j;
  }
  // random shuffling of the ancestors
  return sample(ancestors, ndraws);
}
