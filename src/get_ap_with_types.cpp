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
//' @param status a vector. the statuses on discharge of the patients.
//'
//' @param first_positive_test a vector with the days of the first positive
//' test of the patients.
//'
//' @param source a vector. the sources of the patients.
//'
//' @param nPatients the total number of patients.

//' @export
// [[Rcpp::export]]
NumericVector get_ap_with_types(StringVector status,
                                NumericVector first_positive_test,
                                NumericVector source,
                                int nPatients) {

  NumericVector vector_ap;

  int is_source = 0;

  for(int pt=0; pt<nPatients; pt++) {
    if(status[pt]!="s") {
      if(first_positive_test[pt]==20000) {
        is_source = 0;
        for (int k=0; k<nPatients; k++){
          if (source[k] == pt+1){
            is_source = 1;
          }
        }
        if (is_source == 0 ){
          vector_ap.push_back(pt+1);
        }
      }
    }
  }
  return vector_ap;
}
