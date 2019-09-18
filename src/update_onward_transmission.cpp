#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Get first onward transmission of a patient
//'
//' This function returns the first onward transmission for a given patient.
//'
//'
//' @param pt an integer. the patient for which the first onward transmission
//' is returned.
//'
//' @param source a vector. the sources of the patients.
//'
//' @param colonisation a vector with days of colonisation of the patients.
//'
//' @param n_patients the total number of patients.
//' @export
// [[Rcpp::export]]
int update_onward_transmission(int pt,
                               NumericVector source,
                               NumericVector colonisation,
                               int n_patients){

  int first_onward_transmission = 20000;
  for (int k = 0; k < n_patients; k++){
    if (source[k] == pt){
      if (colonisation[k] < first_onward_transmission){
        first_onward_transmission = colonisation[k];
        }
    }
    }
  return first_onward_transmission;
}
