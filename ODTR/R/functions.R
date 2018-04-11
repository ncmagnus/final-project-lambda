#' Additional Functions
#'
#' Helps odtr to design optimal treatment rule and compute mean outcome
#' @importFrom stats predict glm var qnorm qlogis
#' @importFrom SuperLearner SuperLearner
#' @importFrom hitandrun simplex.sample
#' @export



# SL.correct
# correctly specified QAW regression for simulations (DGP 2 and 3)
SL.correct = function (Y, X, newX, family, obsWeights, model = TRUE, ...) {
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ A + W1 + A:W2 + A:W1:W2, data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}


# dopt.fun
# function that takes as input blip and kappa
# outputs dopt. If kappa is present, computes dopt with resource
# constraints according to kappa = proportion of poeple who can get treatment.
dopt.fun = function(blip, kappa){
  if (is.null(kappa)) {
    dopt = as.numeric(blip > 0)
  } else {
    # estimate S0
    tau = seq(from = min(blip), to = max(blip), by = 0.01) # let tau vary from 0 to 1
    surv = sapply(tau, function(x) mean(blip > x)) #probability that the blip is greater than some varying tau
    # estimate nu
    nu = ifelse(sum(surv<=kappa)==0, 0, min(tau[which(surv <= kappa)])) #the biggest tau such that the survival prob is <= kappa
    # estimate tau0
    tauP = max(c(nu, 0)) # max between nu and 0
    # estimate dopt with RC
    stopifnot(!(blip == tauP & tauP > 0))
    dopt = as.numeric(blip > tauP)
  }
  return(dopt)
}


