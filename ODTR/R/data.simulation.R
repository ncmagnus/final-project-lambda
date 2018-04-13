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

  # Treatment under rule
  if (!is.null(dA) & !is.null(a)){
    stop("Can only have dA or a")
  } else if (is.null(a) & is.null(dA)) {
    A_star = A
  } else if (!is.null(a)){
    A_star = a
  } else if (dA == "simple dynamic") {
    A_star = ifelse(W2 > 0, 1, 0)
  } else if (dA == "ODTR"){
    A_star = as.numeric(blip > 0)
  } else if (dA == "ODTR-RC" & is.null(kappa)){
    stop("If you have dA as ODTR-RC you must specify a kappa")
  } else if (dA == "ODTR-RC"){
    tau = seq(from = min(blip), to = max(blip), by = 0.01) # let tau vary from min blip to max blip
    surv = sapply(tau, function(x) mean(blip > x)) #probability that the blip is greater than some varying tau
    nu = min(tau[which(surv <= kappa)]) #the biggest tau such that the survival prob is <= kappa
    tauP = max(c(nu, 0)) # max between nu and 0
    A_star = as.numeric(blip > tauP)
  }

  # Outcome
  Y_star = as.numeric( U.Y < plogis(2*A_star + 0.7*W1 - 2*A_star*W2 - A_star*W1*W2))

  # Data and target parameter
  O = data.frame(W1, W2, A, A_star, Y, Y_star)

  return(O)


}
