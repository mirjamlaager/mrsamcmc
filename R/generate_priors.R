#' Generate priors.
#'
#' This function generates a list of functions and a list of parameters
#' to be used as priors in \code{run_mcmc()}. If nothing is supplied,
#' default priors are returned.
#'
#' @param f_beta the function to be used as a prior of beta. The default
#' is Cauchy.
#'
#' @param f_b the function to be used as a prior of b. The default
#' is Normal.
#'
#' @param f_s the function to be used as a prior of s. The default
#' is Normal.
#'
#' @param f_rho_1 the function to be used as a prior of rho_1. The default
#' is Uniform.
#'
#' @param f_rho_2 the function to be used as a prior of rho_2. The default
#' is Uniform.
#'
#' @param p_beta a vector with parameters for the prior of beta. The default
#' is c(0, 4).
#'
#' @param p_b a vector with parameters for the prior of b. The default
#' is c(1, 0.5).
#'
#' @param p_s a vector with parameters for the prior of s. The default
#' is c(1, 0.5).
#'
#' @param p_rho_1 a vector with parameters for the prior of rho_1. The default
#' is c(0, 1).
#'
#' @param p_rho_2 a vector with parameters for the prior of rho_2. The default
#' is c(0, 1).
#'
#' @examples
#' ## use default priors
#' generate_priors()
#'
#' ## modify prior for beta
#' generate_priors(f_beta = dnorm, p_beta = c(1, 0.5))


#' @export
generate_priors <- function(f_beta = NULL,
                            f_b = NULL,
                            f_s = NULL,
                            f_rho_1 = NULL,
                            f_rho_2 = NULL,
                            p_beta = NULL,
                            p_b = NULL,
                            p_s = NULL,
                            p_rho_1 = NULL,
                            p_rho_2 = NULL){


  #define default priors
  prior_functions <- list()
  prior_functions$f_beta <- stats::dcauchy
  prior_functions$f_b <- stats::dnorm
  prior_functions$f_s <- stats::dnorm
  prior_functions$f_rho_1 <- stats::dunif
  prior_functions$f_rho_2 <- stats::dunif


  #set default parameters for default priors
  prior_params <- list()
  prior_params$p_beta <- c(0, 4)
  prior_params$p_b <- c(1, 0.5)
  prior_params$p_s <- c(1, 0.5)
  prior_params$p_rho_1 <- c(0, 1)
  prior_params$p_rho_2 <- c(0, 1)


  #modify default prior functions according to user supplied functions

    if (is.null(f_beta) == F){
      if (is.function(f_beta)){
        prior_functions$f_beta <- f_beta
      } else {
          warning(c("The prior supplied for beta is not a function.
                        Using default prior instead."))
        }
    }

  if (is.null(f_b) == F){
    if (is.function(f_b)){
      prior_functions$f_b <- f_b
    } else {
        warning(c("The prior supplied for b is not a function.
                        Using default prior instead."))
      }
  }

  if (is.null(f_s) == F){
    if (is.function(f_s)){
      prior_functions$f_s <- f_s
    } else {
        warning(c("The prior supplied for s is not a function.
                        Using default prior instead."))
      }
  }

  if (is.null(f_rho_1) == F){
    if (is.function(f_rho_1)){
      prior_functions$f_rho_1 <- f_rho_1
    } else {
        warning(c("The prior supplied for rho_1 is not a function.
                        Using default prior instead."))
      }
  }

  if (is.null(f_rho_2) == F){
    if (is.function(f_rho_2)){
      prior_functions$f_rho_2 <- f_rho_2
    } else {
        warning(c("The prior supplied for rho_2 is not a function.
                        Using default prior instead."))
      }
  }




  #modify default parameters according to user supplied values

  if (is.null(p_beta) == F){
    if (is.numeric(p_beta)){
      prior_params$p_beta <- p_beta
    } else {
        warning(c("The parameters supplied for beta are not numeric.
                        Using default values instead."))
      }
  }


  if (is.null(p_b) == F){
    if (is.numeric(p_b)){
      prior_params$p_b <- p_b
    } else {
      warning(c("The parameters supplied for b are not numeric.
                Using default values instead."))
    }
  }

  if (is.null(p_s) == F){
    if (is.numeric(p_s)){
      prior_params$p_s <- p_s
    } else {
      warning(c("The parameters supplied for s are not numeric.
                Using default values instead."))
    }
  }


  if (is.null(p_rho_1) == F){
    if (is.numeric(p_rho_1)){
      prior_params$p_rho_1 <- p_rho_1
    } else {
      warning(c("The parameters supplied for rho_1 are not numeric.
                Using default values instead."))
    }
  }


  if (is.null(p_rho_2) == F){
    if (is.numeric(p_rho_2)){
      prior_params$p_rho_2 <- p_rho_2
    } else {
      warning(c("The parameters supplied for rho_2 are not numeric.
                Using default values instead."))
    }
  }


  priors <- list(prior_functions, prior_params)
  names(priors) <- c("functions", "parameters")

  return(priors)
}
