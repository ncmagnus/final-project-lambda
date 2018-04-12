#' Simulates data to use with odtr or sodtr
#'
#' simulate data involving treatment, covariates and outcome for optimal or
#' simple dynamic treatment
#' @param n size of simulated data
#' @param kappa optional numeric indicating proportion of people who can receive treatment
#' @export
data.simulation = function(n, dA = NULL, a = NULL, kappa = NULL){

  # Unobserved variable
  H = rbinom(n, size = 1, prob = 0.5)

  # Covariates
  W1 = rnorm(n, mean = 0, sd = 1)
  W2 = rnorm(n, mean = 0, sd = 1)
  W3 = rnorm(n, mean = 0, sd = 1)
  W4 = rnorm(n, mean = 0, sd = 1)

  # Treatment
  A = rbinom(n, size = 1, prob = 0.5)

  # Outcome
  Y = rbinom(n, size = 1, prob = 0.5*plogis(1-W1^2 + 3*W2 + 5*W3^2*A - 4.45*A)+0.5*plogis(-0.5- W3 + 2*W1*W2 + 3*abs(W2)*A - 1.5*A))

  # Blip function
  Qbar1 = 0.5*plogis(1-W1^2 + 3*W2 + 5*W3^2*1 - 4.45*1)+0.5*plogis(-0.5- W3 + 2*W1*W2 + 3*abs(W2)*1 - 1.5*1)
  Qbar0 = 0.5*plogis(1-W1^2 + 3*W2 + 5*W3^2*0 - 4.45*0)+0.5*plogis(-0.5- W3 + 2*W1*W2 + 3*abs(W2)*0 - 1.5*0)
  blip = Qbar1 - Qbar0

  # Treatment under rule
  if (!is.null(dA) & !is.null(a)){
    stop("Can only have dA or a")
  } else if (is.null(a) & is.null(dA)) {
    A_star = A
  } else if (!is.null(a)){
    A_star = a
  } else if (dA == "simple dynamic") {
    A_star = ifelse(W1 > 0, 1, 0)
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

  Y_star = rbinom(n, size = 1, prob = 0.5*plogis(1-W1^2 + 3*W2 + 5*W3^2*A_star - 4.45*A_star)+0.5*plogis(-0.5- W3 + 2*W1*W2 + 3*abs(W2)*A_star - 1.5*A_star))

  # Data and target parameter
  O = data.frame(W1, W2, W3, W4, A, A_star, Y, Y_star)

  return(O)

}
