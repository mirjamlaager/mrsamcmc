#include <Rcpp.h>
#include <algorithm>
#include "daily_col_prob_precalculated.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Log likelihood of transmission
//'
//' This function calculates the transmission component of the overall log
//' likelihood.
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
//'
//' @export
// [[Rcpp::export]]
double log_lik_transmission(NumericVector &admissions,
                            NumericVector &colonisations,
                            NumericVector &discharges,
                            StringVector &statuses,
                            NumericMatrix &antibiotics,
                            NumericVector &cp_no_abx,
                            NumericVector &cp_abx,
                            double beta,
                            double b,
                            double s) {

  int n_patients = admissions.size();
  double ll = 0;

  for (int k=0;k<n_patients;k++) {
    if (statuses[k] == "s"){
      //prob avoid infection
      for (int day = admissions[k]; day <= discharges[k]; day++){
        //add prob avoid infection on each day

        ll += daily_col_prob_precalculated(admissions,
                                                   colonisations,
                                                   discharges,
                                                   statuses,
                                                   antibiotics,
                                                   cp_no_abx,
                                                   cp_abx,
                                                   beta,
                                                   b,
                                                   s,
                                                   day,
                                                   k);

      }
    }
    //for those that aquire
    else if (statuses[k] == "a"){
      //if infected on day of admission
      if (colonisations[k] == admissions[k]) {

        ll += log(1.0 - exp(daily_col_prob_precalculated(admissions,
                                                                 colonisations,
                                                                 discharges,
                                                                 statuses,
                                                                 antibiotics,
                                                                 cp_no_abx,
                                                                 cp_abx,
                                                                 beta,
                                                                 b,
                                                                 s,
                                                                 colonisations[k],
                                                                 k)));

      }
      else {
        //prob avoid infection until time interval before infecteed
        for (int day = admissions[k]; day < colonisations[k]; day++){

          ll += daily_col_prob_precalculated(admissions,
                                                     colonisations,
                                                     discharges,
                                                     statuses,
                                                     antibiotics,
                                                     cp_no_abx,
                                                     cp_abx,
                                                     beta,
                                                     b,
                                                     s,
                                                     day,
                                                     k);


        }
        //prob infected on day of infection

        ll+= log(1.0 - exp(daily_col_prob_precalculated(admissions,
                                                                colonisations,
                                                                discharges,
                                                                statuses,
                                                                antibiotics,
                                                                cp_no_abx,
                                                                cp_abx,
                                                                beta,
                                                                b,
                                                                s,
                                                                colonisations[k],
                                                                k)));

      }
    }
  }
  return ll;
}
