#' Compute Simple Dynamic Treatment Rule
#'
#' Designs simple treatment rule and computes mean outcome
#' under this simple rule
#' @param W a dataframe of covariates
#' @param W_for_g a subset of W used to estimate the treatment mechanism
#' @param A a binary vector indicating observed treatment
#' @param a treatment for comparison with A
#' @param Y vector for outcome
#' @param rule vector of treatment (A) under dynamic treatment rule
#' @param QAW.SL.library Super learner library for estimating outcome regression
#' see SuperLearner help file for more information
#' @param estimator.type character indicating estimator type to calculate
#' difference in expected outcome under simple dynamic rule and for comparison group
#' options: gcomp, IPTW, IPTW_DR, TMLE
#' note: selecting IPTW and gcomp preclude confidence intervals
#' default TMLE
#' @importFrom SuperLearner SuperLearner
#' @importFrom stats predict glm qnorm
#' @usage
#' sdtr(W, W_for_g, A, a, Y, rule, QAW.SL.library, estimator.type="TMLE")
#' @export
#

# function that computes gcomp, IPTW, TMLE for simple dynamic txt regime
sdtr = function(W, W_for_g, A, a, Y, rule, QAW.SL.library, estimator.type = "TMLE"){

  #sanity checks

  # check W
  if (!is.data.frame(W)) stop("W should be a dataframe")

  # check W_for_g
  if (!is.data.frame(W_for_g) & !is.vector(W_for_g)) stop("W_for_g should be a dataframe or a vector")

  # check A
  if (!is.vector(A) | length(unique(A)) > 2) stop("A should be a binary vector")

  # check a
  if (!is.vector(a)) stop("a should be a vector for comparison with A")

  # check Y
  if (!(is.vector(Y)) | !is.numeric(Y)) stop("Y should be a numeric vector")

  #check rule
  if (!(is.vector(rule)) | !is.numeric(rule)) stop("rule should be a numeric vector")


  # check QAV.SL.library
  if (!is.character(QAW.SL.library) & !is.list(QAW.SL.library)) stop("QAW.SL.library should be a character vector or a list
                                                                     containing character vectors")

  # check estimator.type
  if (!is.character(estimator.type)) stop("estimator.type should be a character; see help file" )
  if(!estimator.type %in% c("gcomp","TMLE","IPTW","IPTW_DR")){stop("estimator.type must be either 'gcomp', 'TMLE','IPTW', or 'IPTW_DR'")}

  n = length(A)
  family = ifelse(length(unique(Y))>2, "gaussian", "binomial")

  # estimate E[Y|A,W]
  Q.reg = SuperLearner(Y = Y, X = data.frame(A, W), SL.library = QAW.SL.library, family = family)
  # get estimates E[Y|A,W=w]
  QAW = predict(Q.reg, type = "response")$pred
  # get estimates E[Y|A = rule, W=w]
  Qd = predict(Q.reg, newdata = data.frame(W, A = rule), type = "response")$pred
  # get estimates E[Y|A = a, W=w]
  QaW = predict(Q.reg, newdata = data.frame(W, A = a), type = "response")$pred

  ### g-comp ###
  Psi_gcomp = mean(Qd) - mean(QaW)

  # estimate pred. prob. observed exposure, P(A|W)
  g.reg = glm(A ~ ., data = data.frame(A,W_for_g), family = "binomial")
  g1W = predict(g.reg, type = "response")
  gAW = ifelse(A == 1, g1W, 1-g1W)

  ### IPTW ###
  Psi_IPTW = mean(Y*(1/gAW)*as.numeric(A == rule)) - mean(Y*(1/gAW)*as.numeric(A == a))

  ### DR-IPTW ###
  Psi_IPTW_DR = mean(as.numeric(A == rule)/gAW * (Y-Qd) + Qd) - mean(as.numeric(A == a)/gAW * (Y-QaW) + QaW)

  ### TMLE ###
  tmle_objects.d = tmle.fun(A = A, Y = Y, d = rule, Qd = Qd, gAW = gAW, family = family)
  tmle_objects.a = tmle.fun(A = A, Y = Y, d = a, Qd = QaW, gAW = gAW, family = family)
  Psi_TMLE = tmle_objects.d$psi - tmle_objects.a$psi

  # IC IPTW-DR
  IC_IPTW_DR = ((1/gAW)*as.numeric(A == rule)*(Y - Qd) + Qd - mean(as.numeric(A == rule)/gAW * (Y-Qd) + Qd)) - ((1/gAW)*as.numeric(A == a)*(Y - QaW) + QaW - mean(as.numeric(A == a)/gAW * (Y-QaW) + QaW))
  varIC_IPTW_DR = var(IC_IPTW_DR)/n
  CI_IPTW_DR = Psi_IPTW_DR + c(-1,1)*qnorm(0.975)*sqrt(as.numeric(varIC_IPTW_DR))

  # IC TMLE
  varIC_TMLE = var(tmle_objects.d$IC - tmle_objects.a$IC)/n
  CI_TMLE = Psi_TMLE + c(-1,1)*qnorm(0.975)*sqrt(varIC_TMLE)

  if (estimator.type =="gcomp")
  {output <- c(gcomp_simple = Psi_gcomp)}


  if (estimator.type =="IPTW")
  {output <- c(IPTW_simple = Psi_IPTW,
               CI_TMLE_simple = CI_TMLE)}

  if (estimator.type =="IPTW_DR")
  {output <- c(IPTW_simple_DR = Psi_IPTW_DR,
               CI_IPTW_simple_DR = CI_IPTW_DR)}
  else
    {output <- c(TMLE_simple = Psi_TMLE,
                     CI_TMLE_simple = CI_TMLE)}

  return(output)

}


