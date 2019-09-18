#include <Rcpp.h>
#include <algorithm>
#include "col_pop.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Log likelihood of sources
//'
//' This function calculates the source component of the overall log likelihood.
//'
//' @param nPatients the total number of patients.
//'
//' @param nTypes the number of differnt types of colonisations.
//'
//' @param admission a vector with days of admission of the patients.
//'
//' @param colonisation a vector with days of colonisation of the patients.
//'
//' @param source a vector with the sources of the patients.
//'
//' @param discharge a vector with days of discharge of the patients.
//'
//' @param status a vector with the statuses on discharge of the patients.
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @param b the multiplicative effect of antibiotics on transmissibility.
//'
//' @export
// [[Rcpp::export]]
double log_lik_sources(int nPatients,
                 int nTypes,
                 NumericVector admission,
                 NumericVector colonisation,
                 NumericVector source,
                 NumericVector discharge,
                 StringVector status,
                 NumericMatrix antibiotics,
                 double b){

  double ll = 0;
  double C_nabx = 0;
  double C_abx = 0;
  double transmissibility_factor = 1;
  int n_imports = 0;

  for (int pt=0; pt<nPatients; pt++){
    if (status[pt] == "a"){

      NumericVector CP = col_pop(admission,
                                 colonisation,
                                 discharge,
                                 status,
                                 colonisation[pt],
                                 antibiotics);

      int C_nabx = CP(0);
      int C_abx = CP(1);


      transmissibility_factor = 1;
      NumericVector v = antibiotics(source[pt] - 1,_);

      if (std::count(v.begin(),v.end(),colonisation[pt]) > 0){
        transmissibility_factor = b;
        }

      ll = ll + log(transmissibility_factor / (C_nabx + b*C_abx));
    }
    if (status[pt] == "p"){n_imports++;}

  }

  ll = ll + n_imports * log( 1.0 / nTypes);

  return ll;
}
