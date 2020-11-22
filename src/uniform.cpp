#include <Rcpp.h>


double h(double x, double a, double b)
{
  return std::log((1+b*x)/(1+a*x))/(b-a)-x;
}

double hp(double x, double a, double b)
{
  return 1/((1+a*x)*(1+a*x))-1;
}

//' Newton-Raphson algorithm for solving optimal log utiltiy under uniform distributions
//'
//' @param interval the interval for the uniform distribution
//' @param rate the risk-neutral rate
//' @param n the number of iterations to take in the NR method.
//'
//' @description {For finding the root of the function \eqn{h(x)=\log((1+bx)/(1+ax))/(b-a)-x}.}
//' @details {The algorithm is on wikipedia and quite straightforward.}
// [[Rcpp::export]]
double kellyUniform(std::vector<double> interval, double rate = 0.0, unsigned int n = 500)
{
  double a = interval[0];
  double b = interval[1];
  double mu = (a+b)/2.0;
  double v = (b-a)*(b-a)/12.0; // variance btw
  double x0 = (mu-rate)/v;
  std::vector<double> x(n);
  x[0] = x0;
  for(unsigned int i = 1; i < n; i++)
  {
    x[i] = x[i-1]-h(x[i-1], a, b)/hp(x[i-1], a, b);
  }
  return x[n-1];
}
