#' Simulates Data for odtr or sodtr
#'
#' simulate data involving treatment, covariates and outcome for optimal or
#' simple dynamic treatment
#' @param n size of simulated data
#' @param dA optional character indicating how to construct rule:
#' simple dynamic, ODTR or ODTR-RC
#' if ODTR-RC you must include kappa for constraint
#' @param a optional treatment for comparison with A
#' @param kappa optional numeric indicating proportion of people who can receive treatment
#' @importFrom stats rnorm rbinom plogis runif
#' @export
data.simulation = function(n, dA = NULL, a = NULL, kappa = NULL){
  # Unobserved variables
  U.W1 = runif(n, min=0, max=1)
  U.W2 = rnorm(n, mean=0, sd=1)
  U.Y = runif(n, min=0, max=1)

  # Covariates
  W1 = as.numeric( U.W1 < 0.45)
  W2 = 0.75*U.W2

  # Blip function
  Qbar1 = plogis(2*1 + 0.7*W1 - 2*1*W2 - 1*W1*W2)
  Qbar0 = plogis(2*0 + 0.7*W1 - 2*0*W2 - 0*W1*W2)
  blip = Qbar1 - Qbar0

  # Treatment
  A = rbinom(n, size = 1, prob = 0.5)

  # Outcome
  Y = as.numeric(U.Y < plogis(2*A + 0.7*W1 - 2*A*W2 - A*W1*W2))
  # Data and target parameter
  O = data.frame(W1, W2, A, Y)

  return(O)


}
