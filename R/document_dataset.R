#' Simulated patient data.
#'
#' A dataset containing admission and discharge dates, test
#' results and antibiotic use of 450 patients over 90 days.
#'
#'
#' @format a list containing 4 elements:
#' \describe{
#'   \item{patients}{a dataframe with three columns and 450 rows. Each row
#'   corresponds to a patient, the columns are admission, first_positive_test
#'   and discharge in days. If there is no positive test for a patient,
#'   first_positive_test equals 20000.}
#'   \item{test_results_positive}{a matrix with 450 rows and 30 columns. Each
#'   row corresponds to a patient. Entries correspond to the day of a positive
#'   test.}
#'    \item{test_results_negative}{a matrix with 450 rows and 30 columns. Each
#'   row corresponds to a patient. Entries correspond to the day of a negative
#'   test.}
#'    \item{antibiotics}{a matrix with 450 rows and 100 columns. Each
#'   row corresponds to a patient. Entries correspond to the day when an
#'   antibiotic was administered.}
#' }

"patient_data_simulated"
