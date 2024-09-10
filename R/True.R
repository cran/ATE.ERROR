#' True Estimation of Average Treatment Effect
#'
#' This function performs a true estimation of the Average Treatment Effect (ATE) using the generated values for X and Y. 
#' The consistent estimator is calculated as the difference between the expected value of the outcome for the treated group 
#' and the expected value of the outcome for the control group.
#'
#' The expected value for the treated group, E(Y_1), is calculated as the mean of the product of the treatment assignment 
#' and the outcome divided by the estimated propensity score.
#'
#' The expected value for the control group, E(Y_0), is calculated as the mean of the product of the control assignment 
#' and the outcome divided by one minus the estimated propensity score.
#'
#' The propensity score is estimated by applying a logistic regression model to the true values of the covariates and 
#' treatment assignments.
#'
#' @param Y A numeric vector of outcomes.
#' @param A A numeric vector of treatment assignments.
#' @param Z A numeric vector of covariate Z.
#' @param X A numeric vector of covariate X.
#' @return A numeric value representing the estimated treatment effect.
#' @examples
#' library(ATE.ERROR)
#' data(Simulated_data)
#' Y <- Simulated_data$Y
#' A <- Simulated_data$T
#' Z <- Simulated_data$Z
#' X <- Simulated_data$X
#' True_ATE <- True_Estimation(Y, A, Z, X)
#' print(True_ATE)
#' @importFrom stats glm predict
#' @export
True_Estimation <- function(Y, A, Z, X) {
  n <- length(Y)
  logit_model <- glm(A ~ Z + X, family = binomial(link = "logit"))
  e_hat <- predict(logit_model, type = "response")
  E_Y1 <- mean((A * Y) / e_hat)
  E_Y0 <- mean(((1 - A) * Y) / (1 - e_hat))
  tau_hat_true <- E_Y1 - E_Y0
  tau_hat_true
}
