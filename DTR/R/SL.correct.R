#' Parametric QAW Regression
#'
#' Helpful for simulations
#' @param Y vector for outcome varible
#' @param X dataframe for response variables, could be a matrix
#' @param newX dataframe with new response variables for prediction
#' @param family family for glm fit for options see glm
#' @param obsWeights optional weihgts to be applied to glm
#' @param model logic indicating model frame should be included as output default TRUE
#' @param ... additional arguments to be passed to glm
#' @importFrom stats predict glm var qnorm qlogis
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
