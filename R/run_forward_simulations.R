#' Run forward simulations
#'
#'#' @param data a list containing 4 elements:
#' \describe{
#'   \item{patients}{a dataframe with three columns. Each row
#'   corresponds to a patient, the columns are admission, first_positive_test
#'   and discharge. Entries are in days. If there is no positive test for a patient,
#'   first_positive_test equals 20000.}
#'   \item{test_results_positive}{a matrix. Each row corresponds to a patient.
#'   Entries correspond to the day of a positive test.}
#'    \item{test_results_negative}{a matrix. Each row corresponds to a patient.
#'    Entries correspond to the day of a negative test.}
#'    \item{antibiotics}{a matrix. Each row corresponds to a patient.
#'    Entries correspond to the day when an antibiotic was administered.}
#' }
#'
#' @param chains a list with the output from \code{run_mcmc()}. The list contains
#' two dataframes. The first dataframe contains the chains of the parameter
#' values, with iterations stored as specified in the configuration of
#' \code{run_mcmc()}. The second dataframe contains the chains of the statuses,
#' with iterations stored as specified in the configuration of \code{run_mcmc()}.
#'
#' @param n_iter the number of forward simulations.
#'
#' @param seed set a seed for reproducibility of the simulations.
#'
#' @return The function returns a list with three elements:
#' \describe{
#' \item{positive_tests}{the number of positive tests in each simulation.}
#' \item{negative_tests}{the number of negative tests in each simulation.}
#' \item{parameter_samples}{the random samples from the chains that were
#' used for the simulations.}
#' }
#'
#' @export
run_forward_simulations <- function(data, baseline_parameters, intervention_parameters = NULL, n_iter, seed = NULL){

  if (is.null(seed) == F){
    set.seed(seed)
  }

  if (is.null(intervention_parameters)){
    intervention_parameters$modify_transmission = FALSE
    intervention_parameters$decolonisation = FALSE
    intervention_parameters$isolation = FALSE
  }

  # set up patient population
  patients <- data$patients
  n_days <- max(patients$discharge)
  n_patients <- nrow(patients)

  # get days when patients were tested
  days_of_tests <- cbind(data$test_results_positive, data$test_results_negative)
  days_of_tests[which(days_of_tests==0)] <- 20000

  for (row in 1:nrow(days_of_tests)){
    days_of_tests[row,] <- sort(days_of_tests[row,])
    days_of_tests[row, days_of_tests[row,]==20000] <- 0
  }


  #then: generate example: baseline plus three interventions and plot as boxplot -> slide
  #then: complete documentation of run_forward_simulations
  #then: add some plots to introduction vingette (traceplot of parameter)
  #then: write types vignette
  #then: write forward simulation vignette

  # get days when patients were on antibiotics
  days_of_antibiotics <- data$antibiotics

  # get parameter values
  if (intervention_parameters$modify_transmission == F){
    beta <- baseline_parameters$beta
  } else {
    beta <- intervention_parameters$modify_transmission * baseline_parameters$beta
  }
  b <- baseline_parameters$b
  s <- baseline_parameters$s
  rho_1 <- baseline_parameters$rho_1
  rho_2 <- baseline_parameters$rho_2
  phi <- baseline_parameters$phi

  positive_tests_simulated <- rep(0, n_iter)

  for (iter in 1:n_iter){
    positive_tests_simulated_iter <- 0

    patients$colonisation <- 20000
    patients$status <- "s"
    patients$isolation <- 0

    # set up imports
    n_imports <- round(phi*n_patients)
    imported_patients <- sample(1:n_patients,n_imports)
    patients[imported_patients,]$colonisation <- patients[imported_patients,]$admission -1
    patients[imported_patients,]$status <- "p"


    # run transmission, testing and interventions
    for (day in 1:n_days){

      pt_in_ward = which(patients$admission <= day & patients$discharge >= day)

      pt_col <- pt_in_ward[patients[pt_in_ward,]$colonisation < day & patients[pt_in_ward,]$isolation == 0]

      #count patients in ward that are colonised on day and on/off abx
      cp_no_abx <- 0
      cp_abx <- 0

      for (pt in pt_col){
        if (day %in% days_of_antibiotics[pt, ]){cp_abx <- cp_abx + 1
        } else {
          cp_no_abx <- cp_no_abx + 1
        }
      }

      # transmission

      for (pt in pt_in_ward) {
        # patient is susceptible: patient can get infected
        if(patients[pt, ]$status=="s") {

          #check if patient is on abx
          if (day %in% days_of_antibiotics[pt, ]){
            Sabx <- s
          } else {
            Sabx <- 1
          }

          pi_ij <- Sabx * (beta * cp_no_abx + b * beta * cp_abx)
          p_ij <- 1 - exp(-pi_ij)

          if(runif(1) < p_ij) {
            patients[pt, ]$colonisation <- day
            patients[pt, ]$status <- "a"
          }
        }



      }

      # testing
      pt_in_ward_with_positive_test <- 0
      for (pt in pt_in_ward){


        if (day %in% days_of_tests[pt,]){
          #check if patient is positive
          if (patients[pt, ]$colonisation <= day){

            #check if patient is on antibiotics on that day
            if (day %in% days_of_antibiotics[pt, ]){
              #positive patients on abx get a positive result with probability rho_2
              if (runif(1) < rho_2){
                positive_tests_simulated_iter <- positive_tests_simulated_iter + 1
                pt_in_ward_with_positive_test <- append(pt_in_ward_with_positive_test, pt)

              }
            }  else {
              #positive patients not on abx get a positive result with probability rho_1
              if (runif(1) < rho_1){
                positive_tests_simulated_iter <- positive_tests_simulated_iter + 1
                pt_in_ward_with_positive_test <- append(pt_in_ward_with_positive_test, pt)
              }
            }
          }
          }
      }

      pt_in_ward_with_positive_test <- pt_in_ward_with_positive_test[-1]


      # interventions

      # decolonisation
      if (intervention_parameters$decolonisation == TRUE){
        for (pt in pt_in_ward_with_positive_test){
          patients[pt,]$status <- "s"
          patients[pt,]$colonisation <- 20000
        }
      }

      # isolation
      if (intervention_parameters$isolation == TRUE){
        for (pt in pt_in_ward_with_positive_test){
          patients[pt,]$isolation <- 1
        }
      }



    }

  positive_tests_simulated[iter] <- positive_tests_simulated_iter
  }

  output <- positive_tests_simulated
  return(output)
}
