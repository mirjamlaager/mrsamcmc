#include <Rcpp.h>
#include <algorithm>
#include "log_lik_rho.h"
#include "log_lik_transmission_patient_subset.h"
#include "log_lik_importation.h"
#include "log_lik_sources.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Overall log likelihood for one patient
//'
//' This function calculates the summand of the overall log likelihood that
//' changes when the status of one patient is updated. The constant summands
//' are ignored for efficiency.
//'
//'
//' @param admissions a vector with days of admission of the patients.
//'
//' @param discharges a vector with days of discharge of the patients.
//'
//' @param colonisations a vector with days of colonisation of the patients.
//'
//' @param sources a vector with the sources of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param TP_no_abx the number of true positive tests without antibiotics.
//'
//' @param TP_abx the number of true positive tests with antibiotics.
//'
//' @param test_results_negative a matrix where each row corresponds to a
//' patient and the entries correspond to a day with a negative test result.
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
//' @param rho_1 the test sensitivity without antibiotics (baseline).
//'
//' @param rho_2 the test sensitivity with antibiotics.
//'
//' @param nPatients the total number of patients.
//'
//' @param nTypes the number of types of colonisation.
//'
//' @param phi the probability of a patient beeing colonised on admission.
//'
//' @param pt the patient for which the overall log likelihood is computed.
//'
//' @export
// [[Rcpp::export]]
double log_lik_overall_with_types(NumericVector admissions,
                  NumericVector discharges,
                  NumericVector colonisations,
                  NumericVector sources,
                  StringVector statuses,
                  int TP_no_abx,
                  int TP_abx,
                  NumericMatrix test_results_negative,
                  NumericMatrix antibiotics,
                  double beta,
                  double b,
                  double s,
                  double rho_1,
                  double rho_2,
                  int nPatients,
                  int nTypes,
                  double phi,
                  int pt) {

  double ll = log_lik_rho(colonisations,
                              statuses,
                              TP_no_abx,
                              TP_abx,
                              test_results_negative,
                              antibiotics,
                              rho_1,
                              rho_2) +
              log_lik_transmission_patient_subset(admissions,
                                                  colonisations,
                                                  discharges,
                                                  statuses,
                                                  antibiotics,
                                                  beta,
                                                  b,
                                                  s,
                                                  pt) +
                log_lik_sources(nPatients,
                               nTypes,
                               admissions,
                               colonisations,
                               sources,
                               discharges,
                               statuses,
                               antibiotics,
                               b) +
               log_lik_importation(nPatients, statuses, phi);

  return ll;
}

