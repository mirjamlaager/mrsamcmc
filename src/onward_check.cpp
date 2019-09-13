#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Get last possible day for colonisation
//'
//' This function returns the last day on which a patient can become colonised.
//' If an onward transmission depends on that patient, the last possible day is
//' the day before the onward tranmission. If not, it is the minimum of the
//' first positive test and discharge.
//'
//' @param pt the patient for which the onward transmission is calculated.
//'
//' @param admissions a vector with days of admission of the patients.
//'
//' @param colonisations a vector with days of colonisation of the patients.
//'
//' @param discharges a vector with days of discharge of the patients.
//'
//' @param first_positive_test a vector with the days of the first positive
//' test of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param nPatients the total number of patients.
//'
//' @export
// [[Rcpp::export]]
int onward_check(int pt,
                 NumericVector admissions,
                 NumericVector colonisations,
                 NumericVector discharges,
                 NumericVector first_positive_test,
                 StringVector statuses,
                 int nPatients) {

  int maxDate = std::min(discharges[pt - 1], first_positive_test[pt - 1]);

  /* DISABLE FOR NOW
  //for days when infectious (i.e. day after colonised) and an inpatient
  if (colonisation[pt]<maxDate) {
  for (int t=(colonisation[pt]+1); t<maxDate; t++) {
  //get the number of infectious cases
  int n_inf = col_pop_cpp(admission, colonisation, discharge, status, t);
  if(n_inf==1) {
  //if the number is only one - check for aquistion
  int aq_test = 0;
  for(int pt=0; pt<nPatient; pt++) {
  if(status[pt]=="a") {
  if(colonisation[pt]==t) {
  int aq_test = 1;
  break;
  }
  }
  }
  if(aq_test==1) {
  maxDate = (t-1); //i.e. must have been infected in time interval before the onward transmission
  break;
  }
  }
  }
  }
  */
  return maxDate;
}

