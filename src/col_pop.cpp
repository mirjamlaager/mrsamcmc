#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Count the number of colonised patients
//'
//' This function counts the number of colonised patients on and off antibiotics
//' on a given day.
//'
//' @param admissions a vector with days of admission of the patients.
//'
//' @param colonsiations a vector with days of colonisation of the patients.
//'
//' @param discharges a vector with days of discharge of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param day the day for which the colonised number of colonised patients
//' is counted.
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @export
// [[Rcpp::export]]
NumericVector col_pop(NumericVector admissions,
                      NumericVector colonisations,
                      NumericVector discharges,
                      StringVector statuses,
                      int day,
                      NumericMatrix antibiotics) {

  int cp_abx = 0;
  int cp_no_abx = 0;
  int pt = 0;

  NumericVector col_pop;
  NumericVector n_colonised;
  NumericVector v;

  int n_patients = admissions.size();
  // get vector of colonised patients on day
  for (int k=0; k<n_patients; k++) {
    if ((admissions[k] <= day) & (discharges[k] >= day) & (colonisations[k] < day)){
      col_pop.push_back(k);
    }
  }

  // count how many of those are on/off antibiotics on day
  for (int k=0; k<col_pop.size();k++){
    pt = col_pop[k];

    //does abxUse row pt contain the value day?
    v = antibiotics(pt,_);
    if (std::count(v.begin(),v.end(),day) == 0){
      cp_no_abx++;
    }
    else {cp_abx++;}

  }

  n_colonised.push_back(cp_no_abx);
  n_colonised.push_back(cp_abx);
  return n_colonised;
}
