#' Naive Estimation of Average Treatment Effect
#'
#' This function performs a naive estimation of the ATE.
#' This approach gives us the so-called "naive estimate" by ignoring the difference between
#' (X*, Y*) and (X, Y).
#'
#' @param Y_star A numeric vector of outcomes with potential misclassification.
#' @param A A numeric vector of treatment assignments.
#' @param Z A numeric vector of covariate Z.
#' @param X_star A numeric vector of covariate X with measurement error.
#' @return A numeric value representing the estimated treatment effect.
#' @examples
#' library(ATE.ERROR)
#' data(Simulated_data)
#' Y_star <- Simulated_data$Y_star
#' A <- Simulated_data$T
#' Z <- Simulated_data$Z
#' X_star <- Simulated_data$X_star
#' Naive_ATE_XY <- Naive_Estimation(Y_star, A, Z, X_star)
#' print(Naive_ATE_XY)
#' @importFrom stats glm predict
#' @export
Naive_Estimation <- function(Y_star, A, Z, X_star) {
  n <- length(Y_star)
  logit_model <- glm(A ~ Z + X_star, family = binomial(link = "logit"))
  e_hat <- predict(logit_model, type = "response")
  E_Y1 <- mean((A * Y_star) / e_hat)
  E_Y0 <- mean(((1 - A) * Y_star) / (1 - e_hat))
  tau_hat_naive <- E_Y1 - E_Y0
  tau_hat_naive
}
