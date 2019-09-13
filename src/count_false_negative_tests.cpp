#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
//' Count false negative tests
//'
//' This function counts the number of false negative tests conducted with and
//' without antibioitics.
//'
//' @param colonisations a vector with days of colonisation of the patients.
//'
//' @param statuses a vector the statuses on discharge of the patients
//'
//' @param test_results_negative a matrix where each row corresponds to a patient
//' and the entries correspond to a day with a negative test result.
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @export
// [[Rcpp::export]]
NumericVector count_false_negative_tests(NumericVector colonisations,
                                         StringVector statuses,
                                         NumericMatrix test_results_negative,
                                         NumericMatrix antibiotics) {

  int n_patients = test_results_negative.nrow();
  int n_test_results = test_results_negative.ncol();
  int day = 0;
  int FN_no_abx = 0;
  int FN_abx = 0;

  NumericVector v;
  NumericVector FN;

  for (int pt = 0; pt<n_patients; pt++){
    v = antibiotics(pt,_);
    for (int col = 0; col<n_test_results;col++){
      if (test_results_negative(pt,col) != 0){
        day = test_results_negative(pt,col);
        //is the negative test a false negative?
        if (statuses(pt) == "p" | ((statuses(pt) == "a") & (colonisations(pt) <= day))){
          //was the patient on antibiotics on the day of the false negative test?
          if (std::count(v.begin(),v.end(),day) == 0){
            FN_no_abx++;
          }
          else {FN_abx++;}
        }
      }
    }
  }

  FN.push_back(FN_no_abx);
  FN.push_back(FN_abx);
  return(FN);
}
