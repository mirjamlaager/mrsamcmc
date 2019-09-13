#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


//' Get vector of patients eligible for colonisation removal
//'
//' This function returns a vector of patients who are eligible for removal
//' of colonisation. These are all patients who were never tested positive and
//' who are not the only possible source of infection for another patient.
//'
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

//' @export
// [[Rcpp::export]]
NumericVector get_ap(NumericVector admissions,
                     NumericVector colonisations,
                     NumericVector discharges,
                     NumericVector first_positive_test,
                     StringVector statuses,
                     int nPatients) {

  //find days when acquisitions occur
  std::set<int> days_aq;
  for(int pt = 0; pt < nPatients; pt++) {
    if(statuses[pt] == "a") {
      days_aq.insert(colonisations[pt]);
    }
  }

  //get vector of all infectious patients
  NumericVector inf_pt;
  for(int pt = 0; pt < nPatients; pt++) {
    if(statuses[pt] != "s") {
      inf_pt.push_back(pt);
    }
  }

  //find which of these days only one case was infectious
  std::set<int> unique_pt;
  for(int d:days_aq) {
    NumericVector pt_list;
    for(int pt:inf_pt) {
      if( admissions[pt] <= d & discharges[pt] >= d & colonisations[pt] < d) {
        pt_list.push_back(pt);
      }
    }
    if(pt_list.size()==1) {
      unique_pt.insert(pt_list[0]);
    }
  }

  //get final list of patients
  NumericVector vector_ap;
  for (int pt:inf_pt) {
    if(first_positive_test[pt] == 20000) { //get patients without a positive test
      if(std::find(unique_pt.begin(), unique_pt.end(), pt) == unique_pt.end()) { //get patients not in unique_pt
        vector_ap.push_back(pt + 1);
      }
    }
  }
  return vector_ap;
}
