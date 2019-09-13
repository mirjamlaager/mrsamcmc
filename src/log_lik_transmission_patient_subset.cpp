#include <Rcpp.h>
#include <algorithm>
#include "daily_col_prob.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


//' Log likelihood of transmission for one patient
//'
//' This function calculates the summand of the transmission component of the
//' overall log likelihood that changes when the status of one patient is updated.
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
//' @param pt the patient for which the likelihood is computed.

//'
//' @export
// [[Rcpp::export]]
double log_lik_transmission_patient_subset(NumericVector admissions,
                                           NumericVector colonisations,
                                           NumericVector discharges,
                                           StringVector statuses,
                                           NumericMatrix antibiotics,
                                           double beta,
                                           double b,
                                           double s,
                                           int pt) {

  int n_patients = admissions.size();
  double ll = 0;
  int discharge_current_patient = discharges(pt - 1);
  int admission_current_patient = admissions(pt-1);

  for (int k=0;k<n_patients;k++) {
    if ( (admissions[k] <= discharge_current_patient) & (discharges[k] >= admission_current_patient) ){

      // compute contribution of patients who stay susceptible
      if (statuses[k] == "s"){
        //probability of avoiding infection on each day
        for (int day = admissions[k]; day <= discharges[k]; day++){
          ll += daily_col_prob(admissions,
                               colonisations,
                               discharges,
                               statuses,
                               antibiotics,
                               beta,
                               b,
                               s,
                               day,
                               k);
        }
      }

      // compute contribution of patients who acquire
      else if (statuses[k] == "a"){
        // if infected on day of admission
        if (colonisations[k] == admissions[k]) {

          ll += log(1.0 - exp(daily_col_prob(admissions,
                                             colonisations,
                                             discharges,
                                             statuses,
                                             antibiotics,
                                             beta,
                                             b,
                                             s,
                                             colonisations[k],
                                             k)));
        }
        else {
          // probability of avoiding infection on each day until the day before
          // acquisition
          for (int day = admissions[k]; day < colonisations[k]; day++){

            ll += daily_col_prob(admissions,
                                 colonisations,
                                 discharges,
                                 statuses,
                                 antibiotics,
                                 beta,
                                 b,
                                 s,
                                 day,
                                 k);

          }
          // probability of acuqiring on the day of acquisition

          ll+= log(1.0 - exp(daily_col_prob(admissions,
                                            colonisations,
                                            discharges,
                                            statuses,
                                            antibiotics,
                                            beta,
                                            b,
                                            s,
                                            colonisations[k],
                                            k)));
        }
      }
    }
  }
  return ll;
}
