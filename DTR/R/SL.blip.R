#' Computes Blips and Updates SuperLearner and Optimal Rule
#'
#' SuperLearner for estimating the optimal rule
#' @param V subset of W used to construct optimal rule
#' @param W a dataframe of covariates
#' @param A a binary vector indicating observed treatment
#' @param Y vector for outcome
#' @param QAW.reg estimated outcome regression E{Y|A,W}
#' @param gAW estimated treatment mechanism
#' @param QAV.SL.library Super learner library for estimating blip function
#' @param risk.type character of length 1 indicating
#'  risk to be minimized when estimating optimal treatment can be: empirical, TMLE, CV empirical, CV TMLE
#'  default empirical
#' @param family character indicating family for glm
#' options: gaussian or binomial
#' @param grid.size numeric of length 1 indicating size of grid on which to do grid search
#' for weights assigned to algorithms used to estimate optimal treatment
#' default 100
#' @importFrom stats predict glm var qnorm qlogis
#' @importFrom SuperLearner SuperLearner
#' @importFrom hitandrun simplex.sample
#' @export

# SL.blip
# function that takes as input risk type and outputs
# 1. an updated SuperLearner fit of the blip
# 2. the updated coefficients of the blip
# 3. the optimal rule based on the updated SL fit
SL.blip = function(V, W, A, Y, QAW.reg, gAW, QAV.SL.library, risk.type, family, grid.size){

  n = length(A)

  # get estimates E[Y|A = a,W=w]
  QAW = predict(QAW.reg, newdata = data.frame(W, A = A), type = "response")$pred
  # get estimates E[Y|A = 1,W=w]
  Q1W = predict(QAW.reg, newdata = data.frame(W, A = 1), type = "response")$pred
  # get estimates E[Y|A = 0,W=w]
  Q0W = predict(QAW.reg, newdata = data.frame(W, A = 0), type = "response")$pred

  D = (2*A-1)/gAW * (Y-QAW) + Q1W - Q0W # this is same thing as D1 - D0

  #  if(family == "gaussian"){
  # based on mean squared error loss
  SL.blip.init = SuperLearner(Y = D, X = V, SL.library = QAV.SL.library, family = 'gaussian')
  #  } else {
  # based on quasi log-likelihood loss
  #    maxD = boundsY[2] - boundsY[1] - 0.001
  #    minD = boundsY[1] - boundsY[2] + 0.001
  #    nu = function(x) (x-minD)/(maxD-minD)
  #    nu = function(x) (x-min(x))/(max(x)-min(x)) # this transformation is not the same as paper!!!
  #    SL.blip.init = SuperLearner(Y = nu(D), X = V, SL.library = QAV.SL.library, family = "binomial")
  #SL.blip.init = glm(nu(D) ~., data = data.frame(D = D, V), family = "quasi") # TO REPLACE WITH SL
  #  }

  numalgs = length(QAV.SL.library)
  simplex.grid = rbind(SL.blip.init$coef, diag(numalgs), simplex.sample(n = numalgs, N = grid.size)$samples)

  if (risk.type == "empirical") {
    Q.combos = SL.blip.init$Z%*%t(simplex.grid)
    dopt.combos = apply(Q.combos > 0, 2, as.numeric)
    Qdopt.combos = sapply(1:dim(simplex.grid)[1], function(x) predict(QAW.reg, newdata = data.frame(W, A = dopt.combos[,x]), type = "response")$pred)
    risk.combos = sapply(1:dim(simplex.grid)[1], function(x) -mean(((A==dopt.combos[,x])/gAW *(Y-QAW)+Qdopt.combos[,x])))
  } else if (risk.type == "TMLE") {
    Q.combos = SL.blip.init$Z%*%t(simplex.grid)
    dopt.combos = apply(Q.combos > 0, 2, as.numeric)
    Qdopt.combos = sapply(1:dim(simplex.grid)[1], function(x) predict(QAW.reg, newdata = data.frame(W, A = dopt.combos[,x]), type = "response")$pred)
    risk.combos = sapply(1:dim(simplex.grid)[1], function(x) -tmle.fun(A, Y, d = dopt.combos[,x], Qd = Qdopt.combos[,x], gAW = gAW, family = family)$psi)
  } else if (risk.type == "CV empirical"){
    folds = sample(1:10, size = n, replace = T)
    SL.CV.fun = function(i) {
      SL.blip.init.i = SuperLearner(Y = D[folds!=i], X = V[folds!=i,], SL.library = QAV.SL.library, newX = V[folds==i,], family = 'gaussian')
      Q.combo.i = SL.blip.init.i$library.predict%*%t(simplex.grid)
      dopt.combos.i = apply(Q.combo.i > 0, 2, as.numeric)
      Qdopt.combos.i = sapply(1:dim(simplex.grid)[1], function(x) predict(QAW.reg, newdata = data.frame(W[folds==i,], A = dopt.combos.i[,x]), type = "response")$pred)
      risk.combos.i = sapply(1:dim(simplex.grid)[1], function(x) -mean(((A[folds==i]==dopt.combos.i[,x])/gAW[folds==i] *(Y[folds==i]-QAW[folds==i])+Qdopt.combos.i[,x])))
      return(risk.combos.i)
    }
    risk.combos = colMeans(sapply(1:10, SL.CV.fun))
  } else if (risk.type == "CV TMLE"){
    folds = sample(1:10, size = n, replace = T)
    SL.CV.fun = function(i){
      SL.blip.init.i = SuperLearner(Y = D[folds!=i], X = V[folds!=i,], SL.library = QAV.SL.library, newX = V[folds==i,], family = 'gaussian')
      Q.combo.i = SL.blip.init.i$library.predict%*%t(simplex.grid)
      dopt.combos.i = apply(Q.combo.i > 0, 2, as.numeric)
      Qdopt.combos.i = sapply(1:dim(simplex.grid)[1], function(x) predict(QAW.reg, newdata = data.frame(W[folds==i,], A = dopt.combos.i[,x]), type = "response")$pred)
      return(list(dopt.combos.i = dopt.combos.i, Qdopt.combos.i = Qdopt.combos.i))
    }
    SL.CV = lapply(1:10, SL.CV.fun)
    tmle.CV.fun = function(i, x){
      SL.CV.i = SL.CV[[i]]
      dopt.combo.i = SL.CV.i$dopt.combos.i
      Qdopt.combo.i = SL.CV.i$Qdopt.combos.i
      H = (A[folds==i]==dopt.combo.i[,x])/gAW[folds==i]
      return(tmle.fun(A = A[folds==i], Y = Y[folds==i], d = dopt.combo.i[,x], Qd = Qdopt.combo.i[,x], gAW = gAW[folds==i], family = family)$psi)
    }
    risk.combos = colMeans(sapply(1:dim(simplex.grid)[1], function(x) sapply(1:10, function(i) -tmle.CV.fun(i = i, x = x))))
  }

  SL.blip.riskupdate = SL.blip.init

  if (risk.type != "none") {
    SL.blip.riskupdate$alpha = simplex.grid[which.min(risk.combos),]
    SL.blip.riskupdate$SL.predict = SL.blip.init$Z%*%SL.blip.riskupdate$alpha
  }

  return(SL.blip.riskupdate)

}
