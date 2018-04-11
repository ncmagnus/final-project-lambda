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


