#' Generate configuration.
#'
#' This function generates a list of settings to be used as configuration
#' in \code{run_mcmc()}. If nothing is supplied, default settings are returned.
#'
#' @param n_iter total number of iterations including burnin. The
#' default value is 30000.
#'
#' @param n_burnin number of interations in burnin. The default value is 10000
#' or one third of the supplied value for n_iter. During the burnin, the
#' proposal widths are updated to optimise acceptance ratios. The chains
#' returned by run_mcmc include the burnin.
#'
#' @param set_seed logical; if TRUE (default) a seed is set at the beginning
#' of the MCMC iterations to ensure reproducibility of the chain.
#'
#' @param store_params if STORE_ALL (default) all iterations of the chains of
#' the parameters are returned. If a number k is supplied, only every k-th
#' iteration is stored. k = 1 is the same as STORE_ALL.
#'
#' @param store_statuses if STORE_ALL (default) for each patient and each
#' iteration the colonisation time is returned. If a number k is supplied, only
#' every k-th iteration is stored. If AGGREGATED only aggregated (total number of
#' susceptible, acquired, imported) patient data is returned. k = 1 is the
#' same as STORE_ALL.
#'
#' @param infer_covariate_effects if c("b","s","r") (default) covariate effects on
#' transmissibility (b), susceptibility (s) and test sensitivity (r) are inferred.
#' If NONE, no covariate effects are inferred and the values of b, s and r
#' are set to 1. If a vector containing any of the covariate effects parameter
#' b, s and r is supplied, effects on these parameters only are inferred and
#' the remaining parameters are set to 1.
#'
#'
#'
#' @examples
#' ## use default configuration
#' generate_configuration()
#'
#' ## modify number of iterations
#' generate_configuration(n_iter = 1000, n_burnin = 300)
#'
#' ## infer effect of covariates on susceptibility only
#' generate_configuration(infer_covariate_effects = c("s"))
#'


#' @export
generate_configuration <- function(n_iter = NULL,
                                   n_burnin = NULL,
                                   set_seed = NULL,
                                   store_params = NULL,
                                   store_statuses = NULL,
                                   infer_covariate_effects = NULL){


  #define default configuration
  configuration <- list()
  configuration$n_iter <- 30000
  configuration$n_burnin <- 10000
  configuration$set_seed <- TRUE
  configuration$store_params <- "STORE_ALL"
  configuration$store_statuses <- "STORE_ALL"
  configuration$infer_covariate_effects <- c("b", "s", "r")





  #modify default configuration according to user supplied values

  #modify n_iter
  if (is.null(n_iter) == F){
    if (is.numeric(n_iter) == F){
      warning(c("n_iter must be numeric. Using default instead."))
    } else if (n_iter %% 1 != 0){
      warning(c("n_iter must be an integer. Using default instead."))
    } else if (is.null(n_burnin) == F && n_iter < n_burnin){
      warning(c("n_iter must be greater or equal than n_burnin.
                Using default instead."))
    } else {
      configuration$n_iter <- n_iter
      if (is.null(n_burnin) == T){
        configuration$n_burnin <- round(0.3 * n_iter)
        }
    }
  }

  #modify n_burnin
  if (is.null(n_burnin) == F){
    if (is.numeric(n_burnin) == F){
      warning(c("n_burnin must be numeric. Using default instead."))
    } else if (n_burnin %% 1 != 0){
      warning(c("n_burnin must be an integer. Using default instead."))
    } else if (is.null(n_iter) == T & configuration$n_iter < n_burnin){
      warning(c("n_iter must be greater or equal than n_burnin.
                Using default instead."))
    } else if (is.null(n_iter) == F && n_iter < n_burnin){
      warning(c("n_iter must be greater or equal than n_burnin.
                Using default instead."))
    } else {
      configuration$n_burnin <- n_burnin
      }
  }

  #modify set_seed
  if (is.null(set_seed) == F){
    if (is.logical(set_seed) == F){
      warning(c("set_seed must be logical. Using default instead."))
    } else {
      configuration$set_seed <- set_seed
      }
  }

  #modify store_params
  if (is.null(store_params) == F){
    configuration$store_params <- store_params
  }

  #modify store_statuses
  if (is.null(store_statuses) == F){
    configuration$store_statuses <- store_statuses
  }


  #modify infer_covariate_effects
  if (is.null(infer_covariate_effects) == F){
    if (infer_covariate_effects == "NONE"){
      configuration$infer_covariate_effects <- "NONE"
    } else {
      covariate_effects <- c()
      if ("b" %in% infer_covariate_effects){
        append(covariate_effects, "b")
      }

      if ("s" %in% infer_covariate_effects){
        append(covariate_effects, "s")
      }

      if ("r" %in% infer_covariate_effects){
        append(covariate_effects, "r")
      }
      configuration$infer_covariate_effects <- covariate_effects
    }
  }

  return(configuration)
}
