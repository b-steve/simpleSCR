// Template to calculate the negative log-likelihood for a model with
// no random effects. NOT CURRENTLY USED IN THE PACKAGE.
#include <TMB.hpp>
#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_MATRIX(capt);
  DATA_MATRIX(mask_dists);
  DATA_INTEGER(n);
  DATA_INTEGER(n_traps);
  DATA_INTEGER(n_mask);
  DATA_SCALAR(mask_area);
  PARAMETER_VECTOR(pars);
  // Setting a minimum value.
  double dbl_min = 1e-50;
  // Back-transforming parameters.
  Type D = exp(pars(0));
  Type g0 = 1/(1 + exp(-pars(1)));
  Type sigma = exp(pars(2));
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  // Detection probabilities for each mask point.
  vector<Type> prob_det(n_mask);
  // Generating hazard and probability matrices.
  for (int i = 0; i < n_mask; i++){
    for (int j = 0; j < n_traps; j++){
      prob_mat(i, j) = g0*exp(-pow(mask_dists(i, j), 2)/(2*pow(sigma, 2)));
    }
  }
  // The sum of mask probabilities.
  Type sum_prob_det = 0;
  for (int i = 0; i < n_mask; i++){
    Type p_undet = Type(1);
    for (int j = 0; j < n_traps; j++){
      p_undet *= 1 - prob_mat(i, j);
    }
    prob_det(i) = 1 - p_undet;
    sum_prob_det += prob_det(i);
  }
  // Likelihood contributions from capture histories.
  Type log_sum_integrands = 0;
  for (int i = 0; i < n; i++){
    Type integrand = 0;
    for (int j = 0; j < n_mask; j++){
      Type integrand_mask = 0;
      for (int k = 0; k < n_traps; k++){
	integrand_mask += log(pow(prob_mat(j, k), capt(i, k))*pow(1 - prob_mat(j, k), 1 - capt(i, k)) + dbl_min);
      }
      integrand += exp(integrand_mask)*mask_area;
    }
    log_sum_integrands += log(integrand + dbl_min);
  }
  Type f = -log_sum_integrands;
  Type esa = mask_area*sum_prob_det;
  Type n_type = n;
  f -= n*log(D) - D*esa - lgamma(n_type + 1);
  return f;
}
