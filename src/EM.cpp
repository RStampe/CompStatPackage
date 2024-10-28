#include <Rcpp.h>
using namespace Rcpp;

double square(double x) {
  return x * x;
}

struct parameters {
  double mu;
  double sigma2;
  double nu;

  // Constructor
  parameters(double mu_, double sigma2_, double nu_)
    : mu(mu_), sigma2(sigma2_), nu(nu_) {}

  // Easy to return results to r
  List toList() const {
    return List::create(Named("mu") = mu,
                        Named("sigma2") = sigma2,
                        Named("nu") = nu);
  }

  // Easy to get parameters from r
  // static keyword ensure that the function can be called outside of an instance
  static parameters fromList(const List& parList) {
    return parameters(parList["mu"], parList["sigma2"], parList["nu"]);
  }

  double squaredDifference(const parameters& other) const {
    return square(mu - other.mu) +
      square(sigma2 - other.sigma2) +
      square(nu - other.nu);
  }

  double squaredNorm() const {
    return square(mu) + square(sigma2) + square(nu);
  }


};


struct setup {
  NumericVector x;

  setup(NumericVector input_x) : x(input_x) {}

  NumericVector eStep(const parameters& par) {
    double alpha = (par.nu + 1) * par.nu * par.sigma2;
    double beta = par.nu * par.sigma2;

    NumericVector result(x.size());
    for (int i = 0; i < x.size(); i++) {
      result[i] = alpha / (beta + square(x[i] - par.mu));
    }

    return result;
  }

  parameters mStep(NumericVector e, parameters par) {
    double mu_alpha = 0;
    double mu_beta = 0;
    int n = e.size();

    for (int i = 0; i < n; i++) {
      mu_alpha += x[i]*e[i];
    }
    for (int i = 0; i < n; i++) {
      mu_beta += e[i];
    }

    double mu = mu_alpha / mu_beta;

    double sigma2 = 0;
    for (int i = 0; i < n; i++) {
      sigma2 += e[i] * square((x[i] - par.mu))/par.nu;
    }

    sigma2 = sigma2 / static_cast<double>(n);


    return parameters(mu, sigma2, par.nu);
  }
};


bool not_converged(const parameters& current, const parameters& old, double epsilon) {
  double epsilon_squared = epsilon * epsilon;
  double diff_squared = current.squaredDifference(old);
  double threshold = epsilon_squared * std::pow(old.squaredNorm() + epsilon, 2);

  return diff_squared > threshold;
}

//' Perform EM (C++)
//'
//' Finds parameter estimate using EM.
//'
//' @param x Data.
//' @param listOfInitialPar List of initial parameters.
//' @param maxIteration Maximum iterations.
//' @param epsilon Stopping tolerance.
//' @return A list of estimated paramters.
// [[Rcpp::export]]
List perform_em_cpp(NumericVector x, List listOfInitialPar, int maxIteration, double epsilon) {

  parameters initialPar = parameters::fromList(listOfInitialPar);

  setup emSetup = setup(x);

  parameters oldPar = initialPar;

  for(int i = 0; i < maxIteration; i++){
    NumericVector eCal = emSetup.eStep(oldPar);
    parameters newPar = emSetup.mStep(eCal, oldPar);

    if (!not_converged(newPar, oldPar, epsilon)) {
      Rcpp::Rcout << "Converged after " << i + 1 << " iterations.\n";
      break; // Exit the loop if convergence is achieved
    }
    oldPar = newPar;
  }

  return oldPar.toList();
};
