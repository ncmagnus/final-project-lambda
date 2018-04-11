#' Computes TMLE estimate of mean under decision rule
#'
#' Helps odtr to design optimal treatment rule and compute mean outcome
#' @param A a binary vector indicating observed treatment
#' @param Y vector for outcome
#' @param d output from dopt.fun
#' @param Qd Predicted outcome under optimal treatment
#' @param gAW estimate for probability of observed exposure
#' @param family character indicating family for glm
#' options: gaussian or binomial
#' @importFrom stats predict glm var qnorm qlogis
#' @export


# tmle.fun
# function that takes as input A, Y, decision rule, mean under decision rule, gAW, and family
# outputs TMLE estimate of mean under decision rule
tmle.fun = function(A, Y, d, Qd, gAW, family){
  H = (A==d)/gAW
  if(family == "binomial"){
    logit.Qd = qlogis(Qd)
    update = glm(Y~-1+offset(logit.Qd), weights = H, family="quasibinomial")
    Qdopt.star = predict(update, type = "response")
  } else {
    update = glm(Y~-1+offset(Qd), weights = H, family = "gaussian")
    Qdopt.star = predict(update)
  }
  Psi_TMLE = mean(Qdopt.star)
  IC = H*(Y - Qdopt.star) + Qdopt.star - Psi_TMLE
  return(list(psi = Psi_TMLE, IC = IC))
}
