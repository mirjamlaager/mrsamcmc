#' Run forward simulations
#'
#' @param data a list containing 4 elements:
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
run_forward_simulations <- function(data, chains, n_iter, seed = NULL){

  if (is.null(seed) == F){
    set.seed(seed)
  }

  # set up patient population
  patients <- data$patients
  n_days <- max(patients$discharge)
  n_patients <- nrow(patients)

  # get days when patients were tested
  days_of_tests <- cbind(data$test_results_positive, data$test_results_negative)

  # get days when patients were on antibiotics
  days_of_antibiotics <- data$antibiotics

  # sample parameters from mcmc chains
  beta_samples <- sample(chains$parameters$beta, n_iter)
  b_samples <- sample(chains$parameters$b, n_iter)
  s_samples <- sample(chains$parameters$s, n_iter)
  rho_1_samples <- sample(chains$parameters$rho_1, n_iter)
  rho_2_samples <- sample(chains$parameters$rho_2, n_iter)
  phi_samples <- sample(chains$parameters$phi, n_iter)

  positive_tests_simulated <- rep(0, n_iter)
  negative_tests_simulated <- rep(0, n_iter)


  for (iter in 1:n_iter){

  # get parameter values
  beta <- beta_samples[iter]
  b <- b_samples[iter]
  s <- s_samples[iter]
  rho_1 <- rho_1_samples[iter]
  rho_2 <- rho_2_samples[iter]
  phi <- phi_samples[iter]

  # set up imports

  patients$colonisation <- 20000
  patients$status <- "s"

  n_imports <- round(phi*n_patients)
  imported_patients <- sample(1:n_patients,n_imports)
  patients[imported_patients,]$colonisation <- patients[imported_patients,]$admission -1
  patients[imported_patients,]$status <- "p"


  # run infection history
  for (day in 1:n_days){

    pt_in_ward = which(patients$admission <= day & patients$discharge >= day)

    pt_col <- pt_in_ward[patients[pt_in_ward,]$colonisation < day]

    #count patients in ward that are colonised on day and on/off abx
    cp_no_abx <- 0
    cp_abx <- 0

    for (pt in pt_col){
      if (day %in% days_of_antibiotics[pt, ]){cp_abx <- cp_abx + 1
      } else {
        cp_no_abx <- cp_no_abx + 1
        }
    }

    for (pt in pt_in_ward) {
      #patient is susceptible
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

  }

  # run testing
  positive_tests_simulated_iter <- 0
  negative_tests_simulated_iter <- 0

  tests <-  count_tests_forward_simulations(rho_1,rho_2,patients$colonisation,
                                            days_of_antibiotics, days_of_tests,
                                                n_patients)
  positive_tests_simulated_iter <- tests[1]
  negative_tests_simulated_iter <- tests[2]

  # for (pt in 1:n_patients){
  #   print(pt)
  #
  #   #get days of positive and negative test of that patient
  #   days_of_tests_pt <- days_of_tests[pt, days_of_tests[pt,] !=0]
  #
  #   if (length(days_of_tests_pt) > 0){
  #
  #   for (day in days_of_tests){
  #     #check if patient is positive
  #     if (patients[pt, ]$colonisation <= day){
  #
  #       #check if patient is on antibiotics on that day
  #       if (day %in% days_of_antibiotics[pt, ]){
  #         #positive patients on abx get a positive result with probability rho_2
  #         if (runif(1) < rho_2){
  #           positive_tests_simulated_iter <- positive_tests_simulated_iter + 1
  #         } else {
  #             negative_tests_simulated_iter <- negative_tests_simulated_iter + 1
  #             }
  #       }  else {
  #         #positive patients not on abx get a positive result with probability rho_1
  #         if (runif(1) < rho_1){
  #           positive_tests_simulated_iter <- positive_tests_simulated_iter + 1
  #         } else {
  #           negative_tests_simulated_iter <- negative_tests_simulated_iter + 1
  #           }
  #       }
  #     } else {
  #       #negative patients get a negative result with probability 1
  #       negative_tests_simulated_iter <- negative_tests_simulated_iter + 1
  #     }
  #   }
  #   }
  # }

  positive_tests_simulated[iter] <- positive_tests_simulated_iter
  negative_tests_simulated[iter] <- negative_tests_simulated_iter

}

parameter_samples_dataframe <- data.frame(beta_samples,
                                          b_samples,
                                          s_samples,
                                          rho_1_samples,
                                          rho_2_samples,
                                          phi_samples)

output <- list(positive_tests_simulated,
               negative_tests_simulated,
               parameter_samples_dataframe)

names(output) <- c("positive_tests", "negative_tests", "parameter_samples")

return(output)
}
