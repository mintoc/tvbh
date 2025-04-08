// Time-varying Beverton-Holt code for bssm
// uses the template code from Jouni Helske "Non-linear models with bssm" vignette

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// Unknown parameters theta:
// theta(0) = log(H)
// theta(1) = log(R_1)
// theta(2) = log(R_2)

// Function for the prior mean of alpha_1
// [[Rcpp::export]]
arma::vec a1_fn(const arma::vec& theta, const arma::vec& known_params) {
  arma::vec a1(2);
  a1(0) = known_params(0);
  a1(1) = known_params(1);
  return a1;
}
// Function for the prior covariance matrix of alpha_1
// [[Rcpp::export]]
arma::mat P1_fn(const arma::vec& theta, const arma::vec& known_params) {
  arma::mat P1(2, 2, arma::fill::zeros);
  P1(0,0) = known_params(2);
  P1(1,1) = known_params(3);
  return P1;
}

// Function for the observational level standard deviation
// [[Rcpp::export]]
arma::mat H_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat H(1,1);
  H(0, 0) = exp(theta(0));
  return H;
}

// Function for the Cholesky of state level covariance
// [[Rcpp::export]]
arma::mat R0_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat R(2, 2, arma::fill::zeros);
  return R;
}

// R for qvary
// [[Rcpp::export]]
arma::mat R1_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat R(2, 2, arma::fill::zeros);
  R(0, 0) = exp(theta(1));
  return R;
}

// R for pvary
// [[Rcpp::export]]
arma::mat R2_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat R(2, 2, arma::fill::zeros);
  R(1, 1) = exp(theta(1));
  return R;
}


// Z function
// [[Rcpp::export]]
arma::vec Z_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::vec q(1);
  arma::vec p(1);
  arma::vec tau(1);  
  arma::vec fec(1);
  arma::vec tmp(1);
  arma::vec S(1);
  q(0) = exp(alpha(0));
  p(0) = exp(alpha(1));
  tau(0) = known_params(4);
  fec(0) = known_params(5);
  S(0) = known_tv_params(t,0);
  tmp(0) = -log(p(0)/q(0) * (exp(q(0) * tau(0)) - 1.0) + exp(q(0) * tau(0))/(fec(0) * S(0)));
  return tmp;
}

// Jacobian of Z function
// [[Rcpp::export]]
arma::mat Z_gn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::vec lnq(1);
  arma::vec lnp(1);
  arma::vec q(1);
  arma::vec p(1);
  arma::vec tau(1);  
  arma::vec fec(1);
  arma::vec tmp(1);
  arma::vec S(1);
  lnq(0) = alpha(0);
  lnp(0) = alpha(1);  
  q(0) = exp(lnq(0));
  p(0) = exp(lnp(0));
  tau(0) = known_params(4);
  fec(0) = known_params(5);
  S(0) = known_tv_params(t,0);
  //arma::vec ddlnq = arma::vec(1).fill(0.);  
  arma::vec ddlnq(1);
  ddlnq(0) = -(-exp(lnp(0) - lnq(0)) * (exp(exp(lnq(0))) - 1.0) +  exp(lnp(0) + exp(lnq(0))) + (exp(exp(lnq(0)) + lnq(0)))/(S(0) * fec(0)))/(exp(lnp(0) - lnq(0)) * (exp(exp(lnq(0))) - 1.0) + exp(exp(lnq(0)))/(S(0) * fec(0)));
  //arma::vec ddlnp = arma::vec(1).fill(0.);  
  arma::vec ddlnp(1);
  ddlnp(0) = -(exp(lnp(0) - lnq(0)) * (exp(exp(lnq(0))) - 1.0))/(exp(lnp(0) - lnq(0)) * (exp(exp(lnq(0))) - 1.0) + exp(exp(lnq(0)))/(S(0) * fec(0)));
  //
  arma::mat Z_gn(1, 2);
  //
  Z_gn(0, 0) = ddlnq(0);
  Z_gn(0, 1) = ddlnp(0);
  return Z_gn;
}

// T function
// [[Rcpp::export]]
arma::vec T_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::vec alpha_new(2);
  alpha_new(0) = alpha(0);
  alpha_new(1) = alpha(1);
  return alpha_new;
}

// Jacobian of T function
// [[Rcpp::export]]
arma::mat T_gn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, 
  const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat Tg(2, 2);
  Tg(0, 0) = 1.0;
  Tg(0, 1) = 0;
  Tg(1, 0) = 0;
  Tg(1, 1) = 1.0;
  return Tg;
}

// log-prior pdf for theta static
// [[Rcpp::export]]
double log_prior_pdf0(const arma::vec& theta) {
  
  // weakly informative half-N(0, 4) priors. 
  
  // Note that the sampling is on log-scale, 
  // so we need to add jacobians of the corresponding transformations
  // we could also sample on natural scale with check such as
  // if(arma::any(theta < 0)) return -std::numeric_limits<double>::infinity();
  // but this would be less efficient.
  
  // You can use R::dnorm and similar functions, see, e.g.
  // https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  double log_pdf =  
    R::dnorm(exp(theta(0)), 0, 2, 1) +
    arma::accu(theta); //jacobian term
  
  return log_pdf;
}

// log-prior pdf for theta either q or p vary
// [[Rcpp::export]]
double log_prior_pdf12(const arma::vec& theta) {
  
  // weakly informative half-N(0, 4) priors. 
  
  // Note that the sampling is on log-scale, 
  // so we need to add jacobians of the corresponding transformations
  // we could also sample on natural scale with check such as
  // if(arma::any(theta < 0)) return -std::numeric_limits<double>::infinity();
  // but this would be less efficient.
  
  // You can use R::dnorm and similar functions, see, e.g.
  // https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  double log_pdf =  
    R::dnorm(exp(theta(0)), 0, 2, 1) +
    R::dnorm(exp(theta(1)), 0, 1, 1) +
    arma::accu(theta); //jacobian term
  
  return log_pdf;
}

// uniform priors as sometimes the chain wanders into highly negative regions and
// the R session aborts

double log_prior_pdf0_unif(const arma::vec& theta) {
  
  // uniform priors
  
  // Note that the sampling is on log-scale, 
  // so we need to add jacobians of the corresponding transformations
  // we could also sample on natural scale with check such as
  // if(arma::any(theta < 0)) return -std::numeric_limits<double>::infinity();
  // but this would be less efficient.
  
  // You can use R::dnorm and similar functions, see, e.g.
  // https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  double log_pdf =  
    R::dunif(exp(theta(0)), 0.0000001, 2, 1) + 
    arma::accu(theta); //jacobian term
  
  return log_pdf;
}

// log-prior pdf for theta either q or p vary
// [[Rcpp::export]]
double log_prior_pdf12_unif(const arma::vec& theta) {
  
  // weakly informative half-N(0, 4) priors. 
  
  // Note that the sampling is on log-scale, 
  // so we need to add jacobians of the corresponding transformations
  // we could also sample on natural scale with check such as
  // if(arma::any(theta < 0)) return -std::numeric_limits<double>::infinity();
  // but this would be less efficient.
  
  // You can use R::dnorm and similar functions, see, e.g.
  // https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  double log_pdf =  
    R::dunif(exp(theta(0)), 0.0000001, 2, 1) +
    R::dunif(exp(theta(1)), 0.0000001, 1, 1) +
    arma::accu(theta); //jacobian term
  
  return log_pdf;
}



// Create pointers, no need to touch this if
// you don't alter the function names above
// [[Rcpp::export]]
Rcpp::List create_xptrs() {
  
  // typedef for a pointer of nonlinear function of model equation returning vec (T, Z)
  typedef arma::vec (*nvec_fnPtr)(const unsigned int t, const arma::vec& alpha, 
    const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params);
  // typedef for a pointer of nonlinear function returning mat (Tg, Zg, H, R)
  typedef arma::mat (*nmat_fnPtr)(const unsigned int t, const arma::vec& alpha, 
    const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params);
  
  // typedef for a pointer returning a1
  typedef arma::vec (*a1_fnPtr)(const arma::vec& theta, const arma::vec& known_params);
  // typedef for a pointer returning P1
  typedef arma::mat (*P1_fnPtr)(const arma::vec& theta, const arma::vec& known_params);
  // typedef for a pointer of log-prior function
  typedef double (*prior_fnPtr)(const arma::vec& theta);
  
  return Rcpp::List::create(
    Rcpp::Named("a1_fn") = Rcpp::XPtr<a1_fnPtr>(new a1_fnPtr(&a1_fn)),
    Rcpp::Named("P1_fn") = Rcpp::XPtr<P1_fnPtr>(new P1_fnPtr(&P1_fn)),
    Rcpp::Named("Z_fn") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&Z_fn)),
    Rcpp::Named("H_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&H_fn)),
    Rcpp::Named("T_fn") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&T_fn)),
    Rcpp::Named("R0_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R0_fn)),
    Rcpp::Named("R1_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R1_fn)),
    Rcpp::Named("R2_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R2_fn)),
    Rcpp::Named("Z_gn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&Z_gn)),
    Rcpp::Named("T_gn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&T_gn)),
    Rcpp::Named("log_prior_pdf0") = Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf0)),
    Rcpp::Named("log_prior_pdf12") = Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf12)),
    Rcpp::Named("log_prior_pdf0_unif") = Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf0_unif)),
    Rcpp::Named("log_prior_pdf12_unif") = Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf12_unif)));
}
