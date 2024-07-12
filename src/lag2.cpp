#include <RcppEigen.h>
#include <math.h> 
using namespace Rcpp;
using namespace std;
// using namespace arma;

// [[Rcpp::export]]
NumericVector lag2(NumericVector s, int n, double b, String method)
{
  NumericVector w(n);
  std::fill(w.begin(), w.end(), 0.);
  if(method == "bartlett")
  {
    for (int i = 0; i < b; i++)
    {
      w[i] = 1 - s[i]/b;
    }
  }
  else if(method == "tukey")
  {
    for (int i = 0; i < b; i++)
    {
      w[i] = (1 + cos(3.141593 * s[i]/b))/2 ;
    }
  }
  else
  {
    stop("Invalid method. Only bartlett and tukey allowed");
  }

  return(w);
}
