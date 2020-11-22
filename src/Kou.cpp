#include <Rcpp.h>

double pdf_exp(double x, double lambda)
{
  double f = 0.0;
  if (x < 0.0)
  {
    f = 0.0;
  }
  else if (x >= 0.0)
  {
    f = lambda * std::exp(-lambda * x);
  }
  return f;
}

double mgf_exp(double t, double lambda)
{
  double m = 0.0;
  if (t < lambda)
  {
    m = lambda / (lambda - t);
  } else
  {
    Rcpp::stop("MGF for Exp RVs' argument must be less than lambda");
  }
  return m;
}

double pdf_dkou(double x, double prob, double alpha, double beta, double ku, double kd)
{
  double a = 1.0/alpha;
  double b = 1.0/beta;
  double f = 0.0;
  double y;
  if (x >= ku)
  {
    y = x - ku;
    f = prob * pdf_exp(y, a);
  }
  else if (x < kd)
  {
    y = x - kd;
    f = (1.0 - prob)*pdf_exp(-y, b);

  }
  return f;

}

// [[Rcpp::export]]
double mgf_dkou(double t, double prob, double alpha, double beta, double ku, double kd)
{
  double a = 1.0/alpha;
  double b = 1.0/beta;
  double q = 1.0 - prob;
  double m;
  m = q * std::exp(kd)*mgf_exp(-t, b) + prob * std::exp(ku)*mgf_exp(t, a);
  return m;
}
