# function that computes gcomp, IPTW, TMLE for simple dynamic txt regime
simple = function(W, W_for_g, A, a, Y, rule, QAW.SL.library){

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

  return(c(gcomp_simple = Psi_gcomp,
           IPTW_simple = Psi_IPTW,
           IPTW_simple_DR = Psi_IPTW_DR,
           TMLE_simple = Psi_TMLE,
           CI_IPTW_simple_DR = CI_IPTW_DR,
           CI_TMLE_simple = CI_TMLE))

}


