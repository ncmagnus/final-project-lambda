#' Compute Optimal Dynamic Treatment Rule
#'
#' Designs optimal treatment rule and computes mean outcome
#' under this optimal rule
#' @param W a dataframe of covariates
#' @param W_for_g a subset of W used to estimate the treatment mechanism
#' @param A a binary vector indicating observed treatment
#' @param a treatment for comparison with A
#' @param Y vector for outcome
#' @param V subset of W used to construct optimal rule
#' @param QAW.SL.library Super learner library for estimating outcome regression
#' @param QAV.SL.library Super learner library for estimating blip function
#' @param boundsY bounds for outcome variable Y default: c(0,1)
#' @param risk.type character of length 1 indicating
#'  risk to be minimized when estimating optimal treatment can be: empirical, TMLE, CV empirical, CV TMLE
#'  default empirical
#' @param grid.size numeric of length 1 indicating size of grid on which to do grid search
#' for weights assigned to algorithms used to estimate optimal treatment
#' default 100
#' @param kappa (optional) minimum proportion of individuals who can receive treatment (A=1)
#' @importFrom SuperLearner SuperLearner
#' @importFrom hitandrun simplex.sample
#' @importFrom stats predict glm qnorm
#' @usage
#' odtr(W, W_for_g, A, a, Y, V, QAW.SL.library, QAV.SL.library,
#' boundsY = c(0,1), risk.type="empirical", grid.size = 100, kappa = NULL)
#' @export
#



# function that computes gcomp, IPTW, TMLE for optimal dynamic txt regime
odtr = function(W, W_for_g, A, a, Y, V, QAW.SL.library, QAV.SL.library, boundsY = c(0,1), risk.type="empirical", grid.size = 100, kappa = NULL){

  n = length(A)
  family = ifelse(length(unique(Y))>2, "gaussian", "binomial")

  # estimate E[Y|A,W]
  QAW.reg = SuperLearner(Y = Y, X = data.frame(A, W), SL.library = QAW.SL.library, family = family)

  # estimate pred. prob. observed exposure, P(A|W)=g(A|W)
  g.reg = glm(A ~ ., data = data.frame(A,W_for_g), family = "binomial")
  g1W = predict(g.reg, type = "response")
  gAW = ifelse(A == 1, g1W, 1-g1W)

  # get estimate of blip based on risk type (empirical, TMLE, CV TMLE, CV empirical, none)
  SL.blip.fit = SL.blip(V = V, W = W, A = A, Y = Y, QAW.reg = QAW.reg, QAV.SL.library = QAV.SL.library, gAW = gAW, risk.type = risk.type, family = family, grid.size = grid.size)
  # get estimate of optimal rule based on blip estimate
  blip = SL.blip.fit$SL.predict
  # get dopt
  dopt = dopt.fun(blip = blip, kappa = kappa)

  # estimate Qdopt, ie E[Y|A = opt rule,W]
  Qdopt = predict(QAW.reg, newdata = data.frame(W, A = dopt), type = "response")$pred
  # estimate SOC, ie E[Y|A = a, W]
  QaW = predict(QAW.reg, newdata = data.frame(W, A = a), type = "response")$pred

  ### TMLE ###
  tmle_objects.dopt = tmle.fun(A = A, d = dopt, Y = Y, Qd = Qdopt, gAW = gAW, family = family)
  tmle_objects.a = tmle.fun(A = A, d = a, Y = Y, Qd = QaW, gAW = gAW, family = family)
  Psi_TMLE = tmle_objects.dopt$psi - tmle_objects.a$psi

  # IC TMLE
  varIC_TMLE = var(tmle_objects.dopt$IC - tmle_objects.a$IC)/n
  CI_TMLE = Psi_TMLE + c(-1,1)*qnorm(0.975)*sqrt(varIC_TMLE)

  ### CV-TMLE ###
  folds = sample(1:10, size = n, replace = T)
  CV.TMLE = function(i){
    SL.blip.fit.train = SL.blip(V = V[folds!=i,], W = W[folds!=i,], A = A[folds!=i], Y = Y[folds!=i], QAV.SL.library = QAV.SL.library, QAW.reg = QAW.reg, gAW = gAW[folds!=i], risk.type = risk.type, family = family, grid.size = grid.size)
    blip.test = predict(SL.blip.fit.train, newdata = V[folds==i,])$pred
    dopt.test = dopt.fun(blip.test, kappa)
    Qdopt.test = predict(QAW.reg, newdata = data.frame(W[folds == i,], A = dopt.test), type = "response")$pred
    QaW.test = predict(QAW.reg, newdata = data.frame(W[folds == i,], A = a), type = "response")$pred
    tmle_objects.dopt.test = tmle.fun(A = A[folds == i], Y = Y[folds==i], d = dopt.test, Qd = Qdopt.test, gAW = gAW[folds == i], family = family)
    tmle_objects.a.test = tmle.fun(A = A[folds == i], d = a, Y = Y[folds==i], Qd = QaW.test, gAW = gAW[folds == i], family = family)
    Psi_TMLE.test = tmle_objects.dopt.test$psi - tmle_objects.a.test$psi
    var_IC.test = var(tmle_objects.dopt.test$IC - tmle_objects.a.test$IC)
    return(list(Psi_TMLE.test = Psi_TMLE.test, var_IC.test = var_IC.test))
  }
  CV.TMLE.est = lapply(1:10, CV.TMLE)
  Psi_CV.TMLE = mean(sapply(1:10, function(i) CV.TMLE.est[[i]]$Psi_TMLE.test))
  var_CV.TMLE = mean(sapply(1:10, function(i) CV.TMLE.est[[i]]$var_IC.test))/n
  CI_CV.TMLE = Psi_CV.TMLE + c(-1,1)*qnorm(0.975)*sqrt(var_CV.TMLE)
  return(c(TMLE_ODTR = Psi_TMLE,
           CV.TMLE_ODTR = Psi_CV.TMLE,
           CI_TMLE_ODTR = CI_TMLE,
           CI_CV.TMLE_ODTR = CI_CV.TMLE))

}










