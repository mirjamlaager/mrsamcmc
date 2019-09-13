#include <Rcpp.h>
#include <algorithm>
#include "count_false_negative_tests.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Log likelihood of test sensitivity
//'
//' This function calculates the test sensitivity component of the overall log
//' likelihood.
//'
//' @param colonisations a vector with days of colonisation of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param TP_no_abx the number of true positive tests without antibiotics.
//'
//' @param TP_abx the number of true positive tests with antibiotics.
//'
//' @param test_results_negative a matrix where each row corresponds to a patient
//' and the entries correspond to a day with a negative test result.
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @param rho_1 the test sensitivity without antibiotics (baseline).
//'
//' @param rho_2 the test sensitivity with antibiotics.
//'
//'
//' @export
// [[Rcpp::export]]

double log_lik_rho(NumericVector colonisations,
                   StringVector statuses,
                   int TP_no_abx,
                   int TP_abx,
                   NumericMatrix test_results_negative,
                   NumericMatrix antibiotics,
                   double rho_1,
                   double rho_2) {

  NumericVector false_negatives = count_false_negative_tests(colonisations,
                                                             statuses,
                                                             test_results_negative,
                                                             antibiotics);

  int FN_no_abx = false_negatives[0];
  int FN_abx = false_negatives[1];

  double ll = TP_no_abx * log(rho_1) + FN_no_abx * log(1.0 - rho_1) +
              TP_abx * log(rho_2) + FN_abx*log(1.0 - rho_2);

  return(ll);
}


