#include <Rcpp.h>
#include <algorithm>
#include "col_pop.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Daily probability of acquisition
//'
//' This function calculates probability of acquisition for a given patient
//' on a given day.
//'
//' @param admissions a vector with days of admission of the patients.
//'
//' @param colonisations a vector with days of colonisation of the patients.
//'
//' @param discharges a vector with days of discharge of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @param beta the transmission rate.
//'
//' @param b the multiplicative effect of antibiotics on transmissibility.
//'
//' @param s the multiplicative effect of antibiotics on susceptibility.
//'
//' @param day the day for which the probability of acquisition is computed.
//'
//' @param pt the patient for which the probability of acquisition is computed.

//' @export
// [[Rcpp::export]]
double daily_col_prob(NumericVector admissions,
                      NumericVector colonisations,
                      NumericVector discharges,
                      StringVector statuses,
                      NumericMatrix antibiotics,
                      double beta,
                      double b,
                      double s,
                      int day,
                      int pt) {


  NumericVector CP = col_pop(admissions,
                             colonisations,
                             discharges,
                             statuses,
                             day,
                             antibiotics);

  int C_no_abx = CP[0];
  int C_abx = CP[1];

  double abx_factor = 1;
  NumericVector v = antibiotics(pt,_);

  if (std::count(v.begin(),v.end(),day) > 0){abx_factor = s;}

  double p_ij = abx_factor * (-beta * C_no_abx - b * beta * C_abx);

  return p_ij;
}

