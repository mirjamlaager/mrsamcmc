#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Compute the probability of acquisition
//'
//' This function calculates the probability of acquisition for a given day
//' and patient using precalculated numbers of colonised patients on and off
//' antibiotics on the given day.
//'
//' @param admissions a vector with days of admission of the patients.
//'
//' @param colonsiations a vector with days of colonisation of the patients.
//'
//' @param discharges a vector with days of discharge of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @param cp_no_abx the number of colonised patients not on antibiotics
//' on the given day.
//'
//' @param cp_abx the number of colonised patients on antibiotics
//' on the given day.
//'
//' @param beta the transmission rate.
//'
//' @param b the multiplicative effect of antibiotics on transmissibility.
//'
//' @param s the multiplicative effect of antibiotics on susceptibility.
//'
//' @param day the day for which the probability of acquisition is calculated.
//'
//' @param pt the patient for which the probability of acquisition is calculated.
//'
//'
//' @export
// [[Rcpp::export]]
double daily_col_prob_precalculated(NumericVector admissions,
                                    NumericVector colonisations,
                                    NumericVector discharges,
                                    StringVector statuses,
                                    NumericMatrix antibiotics,
                                    NumericVector cp_no_abx,
                                    NumericVector cp_abx,
                                    double beta,
                                    double b,
                                    double s,
                                    int day,
                                    int pt) {


  int C_no_abx = cp_no_abx(day - 1);
  int C_abx = cp_abx(day - 1);

  double abx_factor = 1;
  NumericVector v = antibiotics(pt,_);

  if (std::count(v.begin(), v.end(),day) > 0){abx_factor = s;}

  double p_ij = abx_factor * (-beta * C_no_abx - b * beta * C_abx);

  return p_ij;
}

