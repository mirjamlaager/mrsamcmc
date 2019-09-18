#' Run data augmented MCMC
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
#' @param priors a list containing two lists, functions and params, defining
#' the functions and parameters used as priors. If no input is supplied,
#' default priors are used. Default priors are described in the documentation
#' of \code{generate_priors()}. Custom priors can be generated using
#' \code{generate_priors()}.
#'
#'
#' @param configuration a list of settings. If no input is supplied, default
#' settings are used. Default settings are described in the documentation of
#' \code{generate_configuration()}. Custom settings can be generated using
#' \code{generate_configuration()}.
#
#' @return The function returns a list of two dataframes. The first dataframe
#' contains the chains of the parameter values, with iterations stored as
#' specified in the configuration. The second dataframe contains the chains
#' of the statuses, with iterations stored as specified in the configuration.
#'
#'
#' @examples
#' \dontrun{
#' #run example dataset with default priors and default settings.
#' run_mcmc()
#'
#' #run example dataset with custom priors and custom settings.
#' my_priors <- generate_priors(f_beta = dnorm, p_beta = c(0.02,0.001))
#' my_configuration <- generate_configuration(n_iter = 300, n_burnin = 100)
#' run_mcmc(data = patient_data_simulated, priors = my_priors,
#' configuration = my_configuration)
#'}
#' @export
run_mcmc <- function(data = NULL, priors = NULL, configuration = NULL){

  if (is.null(data) == T){
    data <- patient_data_simulated
  }

  if (is.null(priors) == T){
    priors <- generate_priors()
  }

  if (is.null(configuration) == T){
    configuration <- generate_configuration()
  }

  if (configuration$set_seed == T){
    set.seed(42)
  }

  # initialise chain

  patients <- data$patients
  test_results_positive <- data$test_results_positive
  test_results_negative <- data$test_results_negative
  antibiotics <- data$antibiotics

  n_patients <- nrow(patients)
  n_days <- max(patients$discharge)

  n_iter <- configuration$n_iter
  n_burnin <- configuration$n_burnin

  if ("b" %in% configuration$infer_covariate_effects){
    update_b <- TRUE
  } else {
    update_b <- FALSE
  }

  if ("s" %in% configuration$infer_covariate_effects){
    update_s <- TRUE
  } else {
    update_s <- FALSE
  }

  if ("r" %in% configuration$infer_covariate_effects){
    update_r <- TRUE
  } else {
    update_r <- FALSE
  }


  # initialise augmented data columns.
  # status: s = never infected, a = acquired, p=imported
  patients$colonisation <- 20000
  patients$status <- "s"

  # set up imports: initially as all those infected on the first day
  index <- which(patients$first_positive_test != 20000 &
                 patients$first_positive_test == patients$admission)

  if (length(index) > 0){
    patients[index, ]$colonisation <- patients[index, ]$admission - 1
    patients[index, ]$status <- "p"
  }

  # set up aquisitions: initially as any infected after the first day
  index <- which(patients$first_positive_test != 20000 &
                 patients$first_positive_test > patients$admission)

  if (length(index) > 0){
    patients[index, ]$colonisation <- patients[index, ]$first_positive_test
    patients[index, ]$status <- "a"
  }

  # set up patients who acquire on a day where there was no other
  # colonised patient as imported
  index <- which(patients$status == "a")
  for (pt in index){
    col_time <- patients[pt, ]$colonisation
    pt_in_ward_at_col_time <- which(patients$admission <= col_time &
                                      patients$discharge >= col_time)

    if (sum(patients[pt_in_ward_at_col_time, ]$status == "p") == 0){
      patients[pt, ]$status <- "p"
      patients[pt, ]$colonisation <- patients[pt, ]$admission - 1
    }
  }


  # initial values of the chain
  beta_cur <- 0.1 * runif(1)
  b_cur <- 1
  s_cur <- 1

  rho_1_cur <- 0.5 + 0.5 * runif(1)
  rho_1_logit_cur <- logit(rho_1_cur)

  rho_2_cur <- rho_1_cur
  rho_2_logit_cur <- logit(rho_2_cur)

  phi_cur <- 0.5 * runif(1)

  # proposal densities. The proposal distributions are normal distributions
  # with mean = current value and variance = sigma.
  # The variance of the proposal distributions are adapted during the burnin
  # to optimise acceptance rates.

  sigma_beta <- 0.05
  sigma_b <- 0.5
  sigma_s <- 0.5
  sigma_rho_1 <- 0.05
  sigma_rho_2 <- 0.05

  # prior functions
  prior_function_beta <- priors$functions$f_beta
  prior_function_b <- priors$functions$f_b
  prior_function_s <- priors$functions$f_s
  prior_function_rho_1 <- priors$functions$f_rho_1
  prior_function_rho_2 <- priors$functions$f_rho_2


  # prior parameters
  prior_param_1_beta <- priors$parameters$p_beta[1];
  prior_param_2_beta <- priors$parameters$p_beta[2];
  prior_param_1_b <- priors$parameters$p_b[1];
  prior_param_2_b <- priors$parameters$p_b[2];
  prior_param_1_s <- priors$parameters$p_s[1];
  prior_param_2_s <- priors$parameters$p_s[2];
  prior_param_1_rho_1 <- priors$parameters$p_rho_1[1];
  prior_param_2_rho_1 <- priors$parameters$p_rho_1[2];
  prior_param_1_rho_2 <- priors$parameters$p_rho_2[1];
  prior_param_2_rho_2 <- priors$parameters$p_rho_2[2];




  # initialise accepted moves counting for adaptive mcmc
  N_accepted_beta <- 0
  N_accepted_b <- 0
  N_accepted_s <- 0
  N_accepted_rho_1 <- 0
  N_accepted_rho_2 <- 0

  # tuning parameter for the status updates
  # w is the probability of proposing an importation (as opposed to proposing
  # an acquisistion). The choice of w should not influence convergence, but
  # it does affect the mixing. We use 0.3 based on C.Worby et al.,
  # The Annals of Applied Statistics, 2016.
  # N_local_status_updates is the proportion of patients that get updated
  # in each iteration of the MCMC
  w <- 0.3
  N_local_status_updates <- round(0.01 * n_patients)

  # set up matrices for output storage

  # parameters
  # store each step of the chain
  if ( configuration$store_params == "STORE_ALL" ) {
    output_matrix_params <- matrix(0, nrow = n_iter, ncol = 6)
    store_every_params <- 1
  }

  # store every k-th iteration
  if ( configuration$store_params != "STORE_ALL" ){
    n_store <- round( n_iter / configuration$store_params)
    output_matrix_params <- matrix(0, nrow = n_store, ncol = 6)
    store_every_params <- configuration$store_params
  }

  row_params <- 1

  # statuses
  # store each step of the chain
  if ( configuration$store_statuses == "STORE_ALL" ) {
    output_matrix_statuses <- matrix(0, nrow = n_iter, ncol = n_patients)
    store_every_statuses <- 1
  }

  # store every k-th iteration
  if ( configuration$store_statuses != "STORE_ALL" &&
       configuration$store_statuses != "AGGREGATED" ){
    n_store <- round( n_iter / configuration$store_statuses)
    output_matrix_statuses <- matrix(0, nrow = n_store, ncol = n_patients)
    store_every_statuses <- configuration$store_statuses
  }

  # store aggregated values only
  if (configuration$store_statuses == "AGGREGATED"){
    output_matrix_statuses <- matrix(0, nrow = n_iter, ncol = 3)
    store_every_statuses <- 1
  }

  row_statuses <- 1


  # precalculate the number of true positive tests.
  # since we assume perfect specificity, these numbers are not going to change
  TP <- count_true_positive_tests(test_results_positive, antibiotics)
  TP_no_abx <- TP[1]
  TP_abx <- TP[2]

  CP_no_abx <- rep(0, n_days)
  CP_abx <- rep(0, n_days)

  ### run mcmc ###

   for (iter in 1:n_iter){



    # Adaptive MCMC. During the burnin we update the parameter tuning to
    # optimise the acceptance rates. We aim for an acceptance rate
    # between 0.1 and 0.3. Optimal acceptance rate would be 0.234, as proposed
    # in Roberts, G. O., A. Gelman, and W. R. Gilks. Ann. Appl. Probab. 1997.
    if (iter > 100 & iter < n_burnin & iter %% 100 == 0){

      # adapt proposal width for beta
      if (N_accepted_beta < 10){
        sigma_beta <- 0.8 * sigma_beta;
      }
      else if (N_accepted_beta > 30){
        sigma_beta <- 1.25 * sigma_beta;
      }
      N_accepted_beta <- 0;

      # adapt proposal width for b
      if (N_accepted_b < 10){
        sigma_b <- 0.8 * sigma_b;
      }
      else if (N_accepted_b > 30){
        sigma_b <- 1.25 * sigma_b;
      }
      N_accepted_b <- 0;

      # adapt proposal width for s
      if (N_accepted_s < 10){
        sigma_s <- 0.8 * sigma_s;
      }
      else if (N_accepted_s > 30){
        sigma_s <- 1.25 * sigma_s;
      }
      N_accepted_s <- 0;

      # adapt proposal width for rho_1
      if (N_accepted_rho_1 < 10){
        sigma_rho_1 <- 0.8 * sigma_rho_1;
      }
      else if (N_accepted_rho_1 > 30){
        sigma_rho_1 <- 1.25 * sigma_rho_1;
      }
      N_accepted_rho_1 <- 0;

      # adapt proposal width for rho_2
      if (N_accepted_rho_2 < 10){
        sigma_rho_2 <- 0.8 * sigma_rho_2;
      }
      else if (N_accepted_rho_2 > 30){
        sigma_rho_2 <- 1.25 * sigma_rho_2;
      }
      N_accepted_rho_2 <- 0;
    }


    # precalculate colonised population. The number of colonsied patients only
    # changes in the status updates part of the MCMC. To improve performance
    # these numbers can be precalculated and used for all parameter updates.
    for (day in 1:n_days){
      CP <- col_pop(patients$admission,
                    patients$colonisation,
                    patients$discharge,
                    patients$status,
                    day,
                    antibiotics)
      CP_no_abx[day] <- CP[1]
      CP_abx[day] <- CP[2]
    }


    log_lik_cur <- log_lik_transmission(patients$admission,
                                        patients$colonisation,
                                        patients$discharge,
                                        patients$status,
                                        antibiotics,
                                        CP_no_abx,
                                        CP_abx,
                                        beta_cur,
                                        b_cur,
                                        s_cur)


    # update beta

    beta_cand <- rnorm(1, beta_cur, sigma_beta)

    if (beta_cand > 0){
      log_lik_cand <- log_lik_transmission(patients$admission,
                                           patients$colonisation,
                                           patients$discharge,
                                           patients$status,
                                           antibiotics,
                                           CP_no_abx,
                                           CP_abx,
                                           beta_cand,
                                           b_cur,
                                           s_cur)

      if (log_lik_cand > - Inf) {
        log_prior_cur <- prior_function_beta(beta_cur,
                                             prior_param_1_beta,
                                             prior_param_2_beta,
                                             TRUE)

        log_prior_cand <- prior_function_beta(beta_cand,
                                              prior_param_1_beta,
                                              prior_param_2_beta,
                                              TRUE)

        log_proposal_ratio <- log_lik_cand + log_prior_cand -
          (log_lik_cur + log_prior_cur)

        if (log(runif(1, 0, 1)) < log_proposal_ratio){
          beta_cur <- beta_cand
          N_accepted_beta <- N_accepted_beta + 1
          log_lik_cur <- log_lik_cand
        }
      }
    }



    # update b
    if (update_b == TRUE){
    b_cand <- rnorm(1, b_cur, sigma_b)

    if (b_cand > 0){
      log_lik_cand <- log_lik_transmission(patients$admission,
                                           patients$colonisation,
                                           patients$discharge,
                                           patients$status,
                                           antibiotics,
                                           CP_no_abx,
                                           CP_abx,
                                           beta_cur,
                                           b_cand,
                                           s_cur)

      if (log_lik_cand > - Inf) {

        log_prior_cur <- prior_function_b(b_cur,
                                          prior_param_1_b,
                                          prior_param_2_b,
                                          TRUE)


        log_prior_cand <- prior_function_b(b_cand,
                                           prior_param_1_b,
                                           prior_param_2_b,
                                           TRUE)

        log_proposal_ratio <- log_lik_cand + log_prior_cand -
          (log_lik_cur + log_prior_cur)

        if (log(runif(1, 0, 1)) < log_proposal_ratio){
          b_cur <- b_cand
          N_accepted_b <- N_accepted_b + 1
          log_lik_cur <- log_lik_cand
        }
      }
    }
    }


    # update s
    if (update_s == TRUE){

    s_cand <- rnorm(1, s_cur, sigma_s)

    if (s_cand > 0 ){
      log_lik_cand <- log_lik_transmission(patients$admission,
                                           patients$colonisation,
                                           patients$discharge,
                                           patients$status,
                                           antibiotics,
                                           CP_no_abx,
                                           CP_abx,
                                           beta_cur,
                                           b_cur,
                                           s_cand)

      if (log_lik_cand > - Inf) {
        log_prior_cur <- prior_function_s(s_cur,
                                          prior_param_1_s,
                                          prior_param_2_s,
                                          TRUE)


        log_prior_cand <- prior_function_s(s_cand,
                                           prior_param_1_s,
                                           prior_param_2_s,
                                           TRUE)

        log_proposal_ratio <- log_lik_cand + log_prior_cand -
          (log_lik_cur + log_prior_cur)

        if (log(runif(1, 0, 1)) < log_proposal_ratio){
          s_cur <- s_cand
          N_accepted_s <- N_accepted_s + 1
          log_lik_cur <- log_lik_cand
        }
      }
    }
    }



    # update phi. This is a Gibbs update.
    N_colonised_on_admission <- sum(patients$status == "p")
    N_admissions <- nrow(patients)
    phi_cur <- rbeta(1, N_colonised_on_admission + 1,
                     N_admissions - N_colonised_on_admission + 1)


    # update rho 1. The test sensitivity rho is between 0 and 1, but likely to
    # be close to 1. To ensure good exploration of the space close to one
    # we logit transform the parameters.
    log_lik_rho_cur <- log_lik_rho(patients$colonisation,
                                   patients$status,
                                   TP_no_abx,
                                   TP_abx,
                                   test_results_negative,
                                   antibiotics,
                                   invlogit(rho_1_logit_cur),
                                   invlogit(rho_2_logit_cur))

    rho_1_logit_cand <- rnorm(1, rho_1_logit_cur, sigma_rho_1)

    log_lik_rho_cand <- log_lik_rho(patients$colonisation,
                                    patients$status,
                                    TP_no_abx,
                                    TP_abx,
                                    test_results_negative,
                                    antibiotics,
                                    invlogit(rho_1_logit_cand),
                                    invlogit(rho_2_logit_cur))


    jacobian_cur <- log(invlogit(rho_1_logit_cur)) +
      log(1 - invlogit(rho_1_logit_cur))

    jacobian_cand <- log(invlogit(rho_1_logit_cand)) +
      log(1 - invlogit(rho_1_logit_cand))

    log_proposal_ratio <- log_lik_rho_cand + jacobian_cand -
      (log_lik_rho_cur + jacobian_cur)

    if (log(runif(1, 0, 1)) < log_proposal_ratio){
      rho_1_logit_cur <- rho_1_logit_cand
      rho_1_cur <- invlogit(rho_1_logit_cur)
      N_accepted_rho_1 <- N_accepted_rho_1 + 1
      log_lik_rho_cur <- log_lik_rho_cand
    }



    # update rho_2. The test sensitivity rho is between 0 and 1, but likely to
    # be close to 1. To ensure good exploration of the space close to one
    # we logit transform the parameters.

    if (update_r == TRUE){
    rho_2_logit_cand <- rnorm(1, rho_2_logit_cur, sigma_rho_2)

    log_lik_rho_cand <- log_lik_rho(patients$colonisation,
                                    patients$status,
                                    TP_no_abx,
                                    TP_abx,
                                    test_results_negative,
                                    antibiotics,
                                    invlogit(rho_1_logit_cur),
                                    invlogit(rho_2_logit_cand))


    jacobian_cur <- log(invlogit(rho_2_logit_cur)) +
      log(1 - invlogit(rho_2_logit_cur))

    jacobian_cand <- log(invlogit(rho_2_logit_cand)) +
      log(1 - invlogit(rho_2_logit_cand))

    log_proposal_ratio <- log_lik_rho_cand + jacobian_cand -
      (log_lik_rho_cur + jacobian_cur)

    if (log(runif(1)) < log_proposal_ratio){
      rho_2_logit_cur <- rho_2_logit_cand
      rho_2_cur <- invlogit(rho_2_logit_cur)
      N_accepted_rho_2 <- N_accepted_rho_2 + 1
    }
    } else {
      rho_2_cur <- rho_1_cur
    }


    # update infection times and statuses. In each status update
    # one of three moves is conducted: 1) change a colonisation time, 2) add
    # a colonisation, 3) remove a colonisation.

    for (k in 1:N_local_status_updates){

      # current number of patients in each status
      N_s <- sum(patients$status == "s")
      N_a <- sum(patients$status == "a")
      N_p <- sum(patients$status == "p")

      # sample one of three moves
      type_of_move_sampling <- sample(3, 1)

      # move of type 1: change a colonisation time
      if (type_of_move_sampling == 1){

        # pick a colonised individual
        pt <- sample(which(patients$status != "s"), 1)

        # create copies of the current colonisation and status
        col_cur <- patients$colonisation
        status_cur <- patients$status
        col_cand <- col_cur
        status_cand <- status_cur
        log_hastings_ratio <- NA

        # define window when patient could have aquired. The last possible day
        # of colonisation is the minimum of the discharge and first positive test.
        # To improve mixing, onward transmissions could also be accounted for
        # when defining the window of acquisition. However, it is not obvious,
        # which is faster: rejecting impossible moves (as is now),
        # or checking for onward transmissions.

        first_possible_acquisition_day <- patients[pt, "admission"]

        last_possible_acquisition_day <- onward_check(pt,
                                                      patients$admission,
                                                      patients$colonisation,
                                                      patients$discharge,
                                                      patients$first_positive_test,
                                                      patients$status,
                                                      n_patients)

        dt <- last_possible_acquisition_day - first_possible_acquisition_day + 1

        # propose an importation, with probabilty w
        if (runif(1, 0, 1) < w){
          if (status_cur[pt] == "p") {
            # p -> p (no change)
            log_hastings_ratio <- 0
          } else {
            # a --> p
            col_cand[pt] <- patients[pt, ]$admission - 1
            status_cand[pt] <- "p"
            log_hastings_ratio <- log( (1 - w) / (w * dt))
          }
        } else {
          # propose an acquisition with probability 1 - w
          # sample a new acquisition day at random from all possible days
          status_cand[pt] <- "a"
          if (first_possible_acquisition_day == last_possible_acquisition_day) {
            col_cand[pt] <- first_possible_acquisition_day
          } else {
            col_cand[pt] <- sample(first_possible_acquisition_day:last_possible_acquisition_day, 1)
          }
          if (status_cur[pt] == "p") {
            # p --> a
            log_hastings_ratio <- log( (dt * w) / (1 - w))
          } else {
            # a --> a
            log_hastings_ratio <- 0
          }
        }


        ll_cur <- log_lik_overall(patients$admission,
                                  patients$discharge,
                                  patients$colonisation,
                                  patients$status,
                                  TP_no_abx,
                                  TP_abx,
                                  test_results_negative,
                                  antibiotics,
                                  beta_cur,
                                  b_cur,
                                  s_cur,
                                  rho_1_cur,
                                  rho_2_cur,
                                  n_patients,
                                  phi_cur,
                                  pt)

        ll_cand <- log_lik_overall(patients$admission,
                                   patients$discharge,
                                   col_cand,
                                   status_cand,
                                   TP_no_abx,
                                   TP_abx,
                                   test_results_negative,
                                   antibiotics,
                                   beta_cur,
                                   b_cur,
                                   s_cur,
                                   rho_1_cur,
                                   rho_2_cur,
                                   n_patients,
                                   phi_cur,
                                   pt)


        log_proposal_ratio <- ll_cand - ll_cur + log_hastings_ratio

        if (log(runif(1, 0, 1)) < log_proposal_ratio){
          ll_cur <- ll_cand
          patients$colonisation[pt] <- col_cand[pt]
          patients$status[pt] <- status_cand[pt]
        }

      }

      # move of type 2: add a colonisation
      if (type_of_move_sampling == 2){

        # get all patients who could have a colonisation added. these are all
        # susceptible patients. The number of eligible patients is relevant
        # for the Hastings ratios.
        vector_s <- which(patients$status == "s")
        n_s <- length(vector_s)

        # get all patients who could have a colonisation removed. these are all
        # patients who have never been tested positive. The number of eligible
        # patients is relevant for the Hastings ratios.
        vector_ap <- get_ap(patients$admission,
                            patients$colonisation,
                            patients$discharge,
                            patients$first_positive_test,
                            patients$status,
                            n_patients)
        n_ap <- length(vector_ap)

        # ensure at least one susceptible patient, and sample one
        if (n_s > 0) {
          if (n_s == 1) {
            pt <- vector_s
          } else {
            pt <- sample(vector_s, 1)
          }

          # create copies of the current colonisation and status
          col_cur <- patients$colonisation
          status_cur <- patients$status
          col_cand <- col_cur
          status_cand <- status_cur
          log_hastings_ratio <- NA


          # define window when patient could have aquired. The last possible day
          # of colonisation is discharge, since we assume perfect specificity and
          # never remove colonisations from patients who were tested positive.
          first_possible_acquisition_day <- patients[pt, ]$admission
          last_possible_acquisition_day <- patients[pt, ]$discharge
          dt <- last_possible_acquisition_day - first_possible_acquisition_day + 1

          # propose an importation, with probabilty w
          if (runif(1) < w) {
            # s -> p
            col_cand[pt] <- patients[pt, ]$admission - 1
            status_cand[pt] <- "p"
            log_hastings_ratio <- log(1 / (n_ap + 1)) - log(w / n_s)
          } else {
            #propose an acquisition
            # s -> a
            status_cand[pt] <- "a"
            if (first_possible_acquisition_day == last_possible_acquisition_day) {
              col_cand[pt] <- first_possible_acquisition_day
            } else {
              col_cand[pt] <- sample(first_possible_acquisition_day:last_possible_acquisition_day, 1)
            }
            log_hastings_ratio <- log(1 / (n_ap + 1)) - log( (1 - w) / (n_s * dt))
          }


          ll_cur <- log_lik_overall(patients$admission,
                               patients$discharge,
                               patients$colonisation,
                               patients$status,
                               TP_no_abx,
                               TP_abx,
                               test_results_negative,
                               antibiotics,
                               beta_cur,
                               b_cur,
                               s_cur,
                               rho_1_cur,
                               rho_2_cur,
                               n_patients,
                               phi_cur,
                               pt)


          ll_cand <- log_lik_overall(patients$admission,
                                patients$discharge,
                                col_cand,
                                status_cand,
                                TP_no_abx,
                                TP_abx,
                                test_results_negative,
                                antibiotics,
                                beta_cur,
                                b_cur,
                                s_cur,
                                rho_1_cur,
                                rho_2_cur,
                                n_patients,
                                phi_cur,
                                pt)

          log_proposal_ratio <- ll_cand - ll_cur + log_hastings_ratio

          if (log(runif(1, 0, 1)) < log_proposal_ratio){
            ll_cur <- ll_cand
            patients$colonisation[pt] <- col_cand[pt]
            patients$status[pt] <- status_cand[pt]
          }
        }
      }


      # move of type 3: remove a colonisation
      if (type_of_move_sampling == 3){

        # get all patients who could have a colonisation added. These are all
        # susceptible patients. The number of eligible patients is relevant
        # for the Hastings ratios.
        vector_s <- which(patients$status == "s")
        n_s <- length(vector_s)


        # get all patients who could have a colonisation removed. These are all
        # patients who have never been tested positive. The number of eligible
        # patients is relevant for the Hastings ratios.
        vector_ap <- get_ap(patients$admission,
                            patients$colonisation,
                            patients$discharge,
                            patients$first_positive_test,
                            patients$status,
                            n_patients)
        n_ap <- length(vector_ap)

        # ensure at least one eligible patient, and sample one
        if (n_ap > 0) {
          if (n_ap == 1) {
            pt <- vector_ap
          } else {
            pt <- sample(vector_ap, 1)
          }


          # create copies of the current colinisation and status
          col_cur <- patients$colonisation
          status_cur <- patients$status
          col_cand <- col_cur
          status_cand <- status_cur
          log_hastings_ratio <- NA

          # define window when patient could have aquired. The last possible day
          # of colonisation is the minimum of the discharge and first positive test.
          # To improve mixing, onward transmissions could also be accounted for
          # when defining the window of acquisition. However, it is not obvious,
          # which is faster: rejecting impossible moves (as is now),
          # or checking for onward transmissions.
          first_possible_acquisition_day <- patients[pt, ]$admission
          last_possible_acquisition_day <- onward_check(pt, patients$admission,
                                                        patients$colonisation,
                                                        patients$discharge,
                                                        patients$first_positive_test,
                                                        patients$status,
                                                        n_patients)

          dt <- last_possible_acquisition_day - first_possible_acquisition_day + 1

          status_cand[pt] <- "s"
          col_cand[pt] <- 20000

          if (status_cur[pt] == "p") {
            # p --> s
            log_hastings_ratio <- log(w / (n_s + 1)) - log(1 / n_ap)
          } else if (status_cur[pt] == "a") {
            # a --> s
            log_hastings_ratio <- log( (1 - w) / ( (n_s + 1) * dt)) - log(1 / n_ap)
          }


          ll_cur <- log_lik_overall(patients$admission,
                               patients$discharge,
                               patients$colonisation,
                               patients$status,
                               TP_no_abx,
                               TP_abx,
                               test_results_negative,
                               antibiotics,
                               beta_cur,
                               b_cur,
                               s_cur,
                               rho_1_cur,
                               rho_2_cur,
                               n_patients,
                               phi_cur,
                               pt)

          ll_cand <- log_lik_overall(patients$admission,
                                patients$discharge,
                                col_cand,
                                status_cand,
                                TP_no_abx,
                                TP_abx,
                                test_results_negative,
                                antibiotics,
                                beta_cur,
                                b_cur,
                                s_cur,
                                rho_1_cur,
                                rho_2_cur,
                                n_patients,
                                phi_cur,
                                pt)


          log_proposal_ratio <- ll_cand - ll_cur + log_hastings_ratio

          if (log(runif(1, 0, 1)) < log_proposal_ratio){
            ll_cur <- ll_cand
            patients$colonisation[pt] <- col_cand[pt]
            patients$status[pt] <- status_cand[pt]
          }
        }
      }
    }

    # store values in chain

    # parameters

    if (iter %% store_every_params == 0){
      output_matrix_params[[row_params, 1]] <- beta_cur
      output_matrix_params[[row_params, 2]] <- b_cur
      output_matrix_params[[row_params, 3]] <- s_cur
      output_matrix_params[[row_params, 4]] <- rho_1_cur
      output_matrix_params[[row_params, 5]] <- rho_2_cur
      output_matrix_params[[row_params, 6]] <- phi_cur
      row_params <- row_params + 1
    }

    if (iter %% store_every_statuses == 0){
      if (configuration$store_statuses == "AGGREGATED"){
        output_matrix_statuses[[row_statuses, 1]] <- N_s
        output_matrix_statuses[[row_statuses, 2]] <- N_a
        output_matrix_statuses[[row_statuses, 3]] <- N_p
      } else {
        output_matrix_statuses[row_statuses, ] <- patients$colonisation
      }
      row_statuses <- row_statuses + 1
    }
   }

  output_dataframe_params <- data.frame(output_matrix_params)
  colnames(output_dataframe_params) <- c("beta","b","s","rho_1","rho_2","phi")

  output_dataframe_statuses <- data.frame(output_matrix_statuses)
  if (configuration$store_statuses == "AGGREGATED"){
    colnames(output_dataframe_statuses) <- c("N_never_infected","N_acquired","N_imported")
  }

  output <- list(output_dataframe_params, output_dataframe_statuses)
  names(output) <- c("parameters","statuses")
  return(output)
}
