#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Count true positive tests.
//'
//' This function counts the number of true positive tests conducted with and
//' without antibioitics. We assume perfect specificity, so each positive test
//' is a true positive test.
//'
//' @param test_results_positive a matrix where each row corresponds to a patient
//' and the entries correspond to a day with a positive test result.
//'
//' @param antibiotics a matrix where each row corresponds to a patient
//' and the entries correspond to a day when antibiotics were administered.
//'
//' @export
// [[Rcpp::export]]
NumericVector count_true_positive_tests(NumericMatrix &test_results_positive,
                                        NumericMatrix &antibiotics) {

  int n_patients = test_results_positive.nrow();
  int n_test_results = test_results_positive.ncol();
  int day = 0;
  int TP_no_abx = 0;
  int TP_abx = 0;

  NumericVector v;
  NumericVector TP;

  for (int pt = 0; pt<n_patients; pt++){
    v = antibiotics(pt,_);
    for (int col = 0; col<n_test_results;col++){
      if (test_results_positive(pt,col) != 0){
        day = test_results_positive(pt,col);
        if (std::count(v.begin(),v.end(),day) == 0){
          TP_no_abx++;
        } else {
          TP_abx++;}
      }
    }
  }

  TP.push_back(TP_no_abx);
  TP.push_back(TP_abx);
  return(TP);
}
