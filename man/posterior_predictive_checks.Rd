% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_predictive_checks.R
\name{posterior_predictive_checks}
\alias{posterior_predictive_checks}
\title{Run forward simulations}
\usage{
posterior_predictive_checks(data, chains, n_iter, seed = NULL)
}
\arguments{
\item{data}{a list containing 4 elements:
\describe{
  \item{patients}{a dataframe with three columns. Each row
  corresponds to a patient, the columns are admission, first_positive_test
  and discharge. Entries are in days. If there is no positive test for a patient,
  first_positive_test equals 20000.}
  \item{test_results_positive}{a matrix. Each row corresponds to a patient.
  Entries correspond to the day of a positive test.}
   \item{test_results_negative}{a matrix. Each row corresponds to a patient.
   Entries correspond to the day of a negative test.}
   \item{antibiotics}{a matrix. Each row corresponds to a patient.
   Entries correspond to the day when an antibiotic was administered.}
}}

\item{chains}{a list with the output from \code{run_mcmc()}. The list contains
two dataframes. The first dataframe contains the chains of the parameter
values, with iterations stored as specified in the configuration of
\code{run_mcmc()}. The second dataframe contains the chains of the statuses,
with iterations stored as specified in the configuration of \code{run_mcmc()}.}

\item{n_iter}{the number of forward simulations.}

\item{seed}{set a seed for reproducibility of the simulations.}
}
\value{
The function returns a list with three elements:
\describe{
\item{positive_tests}{the number of positive tests in each simulation.}
\item{negative_tests}{the number of negative tests in each simulation.}
\item{parameter_samples}{the random samples from the chains that were
used for the simulations.}
}
}
\description{
Run forward simulations
}
