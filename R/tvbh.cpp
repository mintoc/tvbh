#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data
  DATA_VECTOR(ssb);
  DATA_VECTOR(y);
  DATA_SCALAR(tau);
  DATA_SCALAR(fec);
  // parameters
  PARAMETER_VECTOR(lnq);
  PARAMETER(lnsdq);
  PARAMETER_VECTOR(lnp);
  PARAMETER(lnsdp);
  PARAMETER(lnsde); // when no information
  // output where the parameters are at
  vector<Type> q = exp(lnq);
  vector<Type> p = exp(lnp);
  vector<Type> a = (p / q) * (exp(q * tau) - 1.0);
  vector<Type> b = exp(q * tau)/fec;
  Type sdq = exp(lnsdq);
  Type sdp = exp(lnsdp);
  Type sde = exp(lnsde);
  int n = ssb.size();
  // negative log-likelihood
  Type nll = 0.0;
  // process
  for(int i = 1; i < n; i++){
    nll -= dnorm(lnq[i], lnq[i-1], sdq, true);
    nll -= dnorm(lnp[i], lnp[i-1], sdp, true);
  }
  // observations
  vector<Type> yhat(n);
  for(int i = 0; i < n; i++){
    yhat[i] = -log(a[i] + b[i]/ssb[i]);
    nll -= dnorm(y[i], yhat[i], sde, true);
  }
  //
  vector<Type> lninvb = log(1.0/b);
  //
  ADREPORT(yhat);
  ADREPORT(lninvb);
  return nll;
}
