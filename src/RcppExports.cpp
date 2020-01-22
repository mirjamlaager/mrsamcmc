// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// col_pop
NumericVector col_pop(NumericVector admissions, NumericVector colonisations, NumericVector discharges, StringVector statuses, int day, NumericMatrix antibiotics);
RcppExport SEXP _mrsamcmc_col_pop(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP statusesSEXP, SEXP daySEXP, SEXP antibioticsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    rcpp_result_gen = Rcpp::wrap(col_pop(admissions, colonisations, discharges, statuses, day, antibiotics));
    return rcpp_result_gen;
END_RCPP
}
// count_false_negative_tests
NumericVector count_false_negative_tests(NumericVector colonisations, StringVector statuses, NumericMatrix test_results_negative, NumericMatrix antibiotics);
RcppExport SEXP _mrsamcmc_count_false_negative_tests(SEXP colonisationsSEXP, SEXP statusesSEXP, SEXP test_results_negativeSEXP, SEXP antibioticsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type test_results_negative(test_results_negativeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    rcpp_result_gen = Rcpp::wrap(count_false_negative_tests(colonisations, statuses, test_results_negative, antibiotics));
    return rcpp_result_gen;
END_RCPP
}
// count_tests_forward_simulations
NumericVector count_tests_forward_simulations(double rho_1, double rho_2, NumericVector colonisations, NumericMatrix antibiotics, NumericMatrix days_of_tests, int n_patients);
RcppExport SEXP _mrsamcmc_count_tests_forward_simulations(SEXP rho_1SEXP, SEXP rho_2SEXP, SEXP colonisationsSEXP, SEXP antibioticsSEXP, SEXP days_of_testsSEXP, SEXP n_patientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho_1(rho_1SEXP);
    Rcpp::traits::input_parameter< double >::type rho_2(rho_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type days_of_tests(days_of_testsSEXP);
    Rcpp::traits::input_parameter< int >::type n_patients(n_patientsSEXP);
    rcpp_result_gen = Rcpp::wrap(count_tests_forward_simulations(rho_1, rho_2, colonisations, antibiotics, days_of_tests, n_patients));
    return rcpp_result_gen;
END_RCPP
}
// count_true_positive_tests
NumericVector count_true_positive_tests(NumericMatrix& test_results_positive, NumericMatrix& antibiotics);
RcppExport SEXP _mrsamcmc_count_true_positive_tests(SEXP test_results_positiveSEXP, SEXP antibioticsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type test_results_positive(test_results_positiveSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type antibiotics(antibioticsSEXP);
    rcpp_result_gen = Rcpp::wrap(count_true_positive_tests(test_results_positive, antibiotics));
    return rcpp_result_gen;
END_RCPP
}
// daily_col_prob
double daily_col_prob(NumericVector admissions, NumericVector colonisations, NumericVector discharges, StringVector statuses, NumericMatrix antibiotics, double beta, double b, double s, int day, int pt);
RcppExport SEXP _mrsamcmc_daily_col_prob(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP statusesSEXP, SEXP antibioticsSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP, SEXP daySEXP, SEXP ptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    rcpp_result_gen = Rcpp::wrap(daily_col_prob(admissions, colonisations, discharges, statuses, antibiotics, beta, b, s, day, pt));
    return rcpp_result_gen;
END_RCPP
}
// daily_col_prob_precalculated
double daily_col_prob_precalculated(NumericVector admissions, NumericVector colonisations, NumericVector discharges, StringVector statuses, NumericMatrix antibiotics, NumericVector cp_no_abx, NumericVector cp_abx, double beta, double b, double s, int day, int pt);
RcppExport SEXP _mrsamcmc_daily_col_prob_precalculated(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP statusesSEXP, SEXP antibioticsSEXP, SEXP cp_no_abxSEXP, SEXP cp_abxSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP, SEXP daySEXP, SEXP ptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cp_no_abx(cp_no_abxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cp_abx(cp_abxSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    rcpp_result_gen = Rcpp::wrap(daily_col_prob_precalculated(admissions, colonisations, discharges, statuses, antibiotics, cp_no_abx, cp_abx, beta, b, s, day, pt));
    return rcpp_result_gen;
END_RCPP
}
// get_ap
NumericVector get_ap(NumericVector admissions, NumericVector colonisations, NumericVector discharges, NumericVector first_positive_test, StringVector statuses, int nPatients);
RcppExport SEXP _mrsamcmc_get_ap(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP first_positive_testSEXP, SEXP statusesSEXP, SEXP nPatientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type first_positive_test(first_positive_testSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type nPatients(nPatientsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ap(admissions, colonisations, discharges, first_positive_test, statuses, nPatients));
    return rcpp_result_gen;
END_RCPP
}
// get_ap_with_types
NumericVector get_ap_with_types(StringVector status, NumericVector first_positive_test, NumericVector source, int nPatients);
RcppExport SEXP _mrsamcmc_get_ap_with_types(SEXP statusSEXP, SEXP first_positive_testSEXP, SEXP sourceSEXP, SEXP nPatientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type first_positive_test(first_positive_testSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type source(sourceSEXP);
    Rcpp::traits::input_parameter< int >::type nPatients(nPatientsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ap_with_types(status, first_positive_test, source, nPatients));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_importation
double log_lik_importation(int nPatient, StringVector statuses, double phi);
RcppExport SEXP _mrsamcmc_log_lik_importation(SEXP nPatientSEXP, SEXP statusesSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nPatient(nPatientSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_importation(nPatient, statuses, phi));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_overall
double log_lik_overall(NumericVector admissions, NumericVector discharges, NumericVector colonisations, StringVector statuses, int TP_no_abx, int TP_abx, NumericMatrix& test_results_negative, NumericMatrix& antibiotics, double beta, double b, double s, double rho_1, double rho_2, int nPatient, double phi, int pt);
RcppExport SEXP _mrsamcmc_log_lik_overall(SEXP admissionsSEXP, SEXP dischargesSEXP, SEXP colonisationsSEXP, SEXP statusesSEXP, SEXP TP_no_abxSEXP, SEXP TP_abxSEXP, SEXP test_results_negativeSEXP, SEXP antibioticsSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP, SEXP rho_1SEXP, SEXP rho_2SEXP, SEXP nPatientSEXP, SEXP phiSEXP, SEXP ptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type TP_no_abx(TP_no_abxSEXP);
    Rcpp::traits::input_parameter< int >::type TP_abx(TP_abxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type test_results_negative(test_results_negativeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type rho_1(rho_1SEXP);
    Rcpp::traits::input_parameter< double >::type rho_2(rho_2SEXP);
    Rcpp::traits::input_parameter< int >::type nPatient(nPatientSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_overall(admissions, discharges, colonisations, statuses, TP_no_abx, TP_abx, test_results_negative, antibiotics, beta, b, s, rho_1, rho_2, nPatient, phi, pt));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_overall_with_types
double log_lik_overall_with_types(NumericVector admissions, NumericVector discharges, NumericVector colonisations, NumericVector sources, StringVector statuses, int TP_no_abx, int TP_abx, NumericMatrix test_results_negative, NumericMatrix antibiotics, double beta, double b, double s, double rho_1, double rho_2, int nPatients, NumericVector types, NumericVector types_names, NumericVector types_freq, int nTypes, double phi, int pt);
RcppExport SEXP _mrsamcmc_log_lik_overall_with_types(SEXP admissionsSEXP, SEXP dischargesSEXP, SEXP colonisationsSEXP, SEXP sourcesSEXP, SEXP statusesSEXP, SEXP TP_no_abxSEXP, SEXP TP_abxSEXP, SEXP test_results_negativeSEXP, SEXP antibioticsSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP, SEXP rho_1SEXP, SEXP rho_2SEXP, SEXP nPatientsSEXP, SEXP typesSEXP, SEXP types_namesSEXP, SEXP types_freqSEXP, SEXP nTypesSEXP, SEXP phiSEXP, SEXP ptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sources(sourcesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type TP_no_abx(TP_no_abxSEXP);
    Rcpp::traits::input_parameter< int >::type TP_abx(TP_abxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type test_results_negative(test_results_negativeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type rho_1(rho_1SEXP);
    Rcpp::traits::input_parameter< double >::type rho_2(rho_2SEXP);
    Rcpp::traits::input_parameter< int >::type nPatients(nPatientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types(typesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types_names(types_namesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types_freq(types_freqSEXP);
    Rcpp::traits::input_parameter< int >::type nTypes(nTypesSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_overall_with_types(admissions, discharges, colonisations, sources, statuses, TP_no_abx, TP_abx, test_results_negative, antibiotics, beta, b, s, rho_1, rho_2, nPatients, types, types_names, types_freq, nTypes, phi, pt));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_rho
double log_lik_rho(NumericVector colonisations, StringVector statuses, int TP_no_abx, int TP_abx, NumericMatrix test_results_negative, NumericMatrix antibiotics, double rho_1, double rho_2);
RcppExport SEXP _mrsamcmc_log_lik_rho(SEXP colonisationsSEXP, SEXP statusesSEXP, SEXP TP_no_abxSEXP, SEXP TP_abxSEXP, SEXP test_results_negativeSEXP, SEXP antibioticsSEXP, SEXP rho_1SEXP, SEXP rho_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type TP_no_abx(TP_no_abxSEXP);
    Rcpp::traits::input_parameter< int >::type TP_abx(TP_abxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type test_results_negative(test_results_negativeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type rho_1(rho_1SEXP);
    Rcpp::traits::input_parameter< double >::type rho_2(rho_2SEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_rho(colonisations, statuses, TP_no_abx, TP_abx, test_results_negative, antibiotics, rho_1, rho_2));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_sources
double log_lik_sources(int nPatients, NumericVector admission, NumericVector colonisation, NumericVector source, NumericVector discharge, StringVector status, NumericVector types, NumericVector types_names, NumericVector types_freq, int nTypes, NumericMatrix antibiotics, double b);
RcppExport SEXP _mrsamcmc_log_lik_sources(SEXP nPatientsSEXP, SEXP admissionSEXP, SEXP colonisationSEXP, SEXP sourceSEXP, SEXP dischargeSEXP, SEXP statusSEXP, SEXP typesSEXP, SEXP types_namesSEXP, SEXP types_freqSEXP, SEXP nTypesSEXP, SEXP antibioticsSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nPatients(nPatientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type admission(admissionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisation(colonisationSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type source(sourceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharge(dischargeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types(typesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types_names(types_namesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type types_freq(types_freqSEXP);
    Rcpp::traits::input_parameter< int >::type nTypes(nTypesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_sources(nPatients, admission, colonisation, source, discharge, status, types, types_names, types_freq, nTypes, antibiotics, b));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_transmission
double log_lik_transmission(NumericVector& admissions, NumericVector& colonisations, NumericVector& discharges, StringVector& statuses, NumericMatrix& antibiotics, NumericVector& cp_no_abx, NumericVector& cp_abx, double beta, double b, double s);
RcppExport SEXP _mrsamcmc_log_lik_transmission(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP statusesSEXP, SEXP antibioticsSEXP, SEXP cp_no_abxSEXP, SEXP cp_abxSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< StringVector& >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type cp_no_abx(cp_no_abxSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type cp_abx(cp_abxSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_transmission(admissions, colonisations, discharges, statuses, antibiotics, cp_no_abx, cp_abx, beta, b, s));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_transmission_patient_subset
double log_lik_transmission_patient_subset(NumericVector admissions, NumericVector colonisations, NumericVector discharges, StringVector statuses, NumericMatrix antibiotics, double beta, double b, double s, int pt);
RcppExport SEXP _mrsamcmc_log_lik_transmission_patient_subset(SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP statusesSEXP, SEXP antibioticsSEXP, SEXP betaSEXP, SEXP bSEXP, SEXP sSEXP, SEXP ptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type antibiotics(antibioticsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_transmission_patient_subset(admissions, colonisations, discharges, statuses, antibiotics, beta, b, s, pt));
    return rcpp_result_gen;
END_RCPP
}
// onward_check
int onward_check(int pt, NumericVector admissions, NumericVector colonisations, NumericVector discharges, NumericVector first_positive_test, StringVector statuses, int nPatients);
RcppExport SEXP _mrsamcmc_onward_check(SEXP ptSEXP, SEXP admissionsSEXP, SEXP colonisationsSEXP, SEXP dischargesSEXP, SEXP first_positive_testSEXP, SEXP statusesSEXP, SEXP nPatientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type admissions(admissionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisations(colonisationsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type discharges(dischargesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type first_positive_test(first_positive_testSEXP);
    Rcpp::traits::input_parameter< StringVector >::type statuses(statusesSEXP);
    Rcpp::traits::input_parameter< int >::type nPatients(nPatientsSEXP);
    rcpp_result_gen = Rcpp::wrap(onward_check(pt, admissions, colonisations, discharges, first_positive_test, statuses, nPatients));
    return rcpp_result_gen;
END_RCPP
}
// update_onward_transmission
int update_onward_transmission(int pt, NumericVector source, NumericVector colonisation, int n_patients);
RcppExport SEXP _mrsamcmc_update_onward_transmission(SEXP ptSEXP, SEXP sourceSEXP, SEXP colonisationSEXP, SEXP n_patientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type source(sourceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type colonisation(colonisationSEXP);
    Rcpp::traits::input_parameter< int >::type n_patients(n_patientsSEXP);
    rcpp_result_gen = Rcpp::wrap(update_onward_transmission(pt, source, colonisation, n_patients));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mrsamcmc_col_pop", (DL_FUNC) &_mrsamcmc_col_pop, 6},
    {"_mrsamcmc_count_false_negative_tests", (DL_FUNC) &_mrsamcmc_count_false_negative_tests, 4},
    {"_mrsamcmc_count_tests_forward_simulations", (DL_FUNC) &_mrsamcmc_count_tests_forward_simulations, 6},
    {"_mrsamcmc_count_true_positive_tests", (DL_FUNC) &_mrsamcmc_count_true_positive_tests, 2},
    {"_mrsamcmc_daily_col_prob", (DL_FUNC) &_mrsamcmc_daily_col_prob, 10},
    {"_mrsamcmc_daily_col_prob_precalculated", (DL_FUNC) &_mrsamcmc_daily_col_prob_precalculated, 12},
    {"_mrsamcmc_get_ap", (DL_FUNC) &_mrsamcmc_get_ap, 6},
    {"_mrsamcmc_get_ap_with_types", (DL_FUNC) &_mrsamcmc_get_ap_with_types, 4},
    {"_mrsamcmc_log_lik_importation", (DL_FUNC) &_mrsamcmc_log_lik_importation, 3},
    {"_mrsamcmc_log_lik_overall", (DL_FUNC) &_mrsamcmc_log_lik_overall, 16},
    {"_mrsamcmc_log_lik_overall_with_types", (DL_FUNC) &_mrsamcmc_log_lik_overall_with_types, 21},
    {"_mrsamcmc_log_lik_rho", (DL_FUNC) &_mrsamcmc_log_lik_rho, 8},
    {"_mrsamcmc_log_lik_sources", (DL_FUNC) &_mrsamcmc_log_lik_sources, 12},
    {"_mrsamcmc_log_lik_transmission", (DL_FUNC) &_mrsamcmc_log_lik_transmission, 10},
    {"_mrsamcmc_log_lik_transmission_patient_subset", (DL_FUNC) &_mrsamcmc_log_lik_transmission_patient_subset, 9},
    {"_mrsamcmc_onward_check", (DL_FUNC) &_mrsamcmc_onward_check, 7},
    {"_mrsamcmc_update_onward_transmission", (DL_FUNC) &_mrsamcmc_update_onward_transmission, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mrsamcmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
