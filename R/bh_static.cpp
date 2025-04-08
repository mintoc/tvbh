#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data
  DATA_VECTOR(ssb);
  DATA_VECTOR(y);
  DATA_SCALAR(tau);
  DATA_SCALAR(fec);
  // parameters
  PARAMETER(lnq);
  PARAMETER(lnp);
  PARAMETER(lnsde); // when no information
  // output where the parameters are at
  Type q = exp(lnq);
  Type p = exp(lnp);
  Type a = (p / q) * (exp(q * tau) - 1.0);
  Type b = exp(q * tau)/fec;
  Type sde = exp(lnsde);
  int n = ssb.size();
  // negative log-likelihood
  Type nll = 0.0;
  // observations
  vector<Type> yhat(n);
  for(int i = 0; i < n; i++){
    yhat[i] = -log(a + b/ssb[i]);
    nll -= dnorm(y[i], yhat[i], sde, true);
  }
  return nll;
}
