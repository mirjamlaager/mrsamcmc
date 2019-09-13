#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Log likelihood of importation
//'
//' This function calculates the importation component of the overall log
//' likelihood.
//'
//' @param nPatient the total number of patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param TP_no_abx the number of true positive tests without antibiotics.
//'
//' @param phi the probability of a patient beeing colonised on admission.
//'
//' @export
// [[Rcpp::export]]

double log_lik_importation(int nPatient,
                           StringVector statuses,
                           double phi) {
  int imports = 0;
  for (int pt=0; pt<nPatient; pt++) {
    if(statuses[pt] == "p") {
      imports ++;
      }
    }
  double ll = imports * log(phi) + (nPatient - imports) * log(1.0 - phi);
  return ll;
}

