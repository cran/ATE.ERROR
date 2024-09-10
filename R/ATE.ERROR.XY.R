#' ATE.ERROR.XY Function for Estimating Average Treatment Effect (ATE) with Measurement Error in X and Misclassification in Y
#'
#' The `ATE.ERROR.XY` function implements a method for estimating the Average Treatment Effect (ATE) that accounts for both measurement error in covariates and misclassification in the binary outcome variable Y.
#' 
#' @param Y_star Numeric vector. The observed binary outcome variable, possibly misclassified.
#' @param A Numeric vector. The treatment indicator (1 if treated, 0 if control).
#' @param Z Numeric vector. A precisely measured covariate vector.
#' @param X_star Numeric vector. A covariate vector subject to measurement error.
#' @param p11 Numeric. The probability of correctly classified Y given Y = 1.
#' @param p10 Numeric. The probability of misclassified Y given Y = 0.
#' @param sigma_epsilon Numeric. The covariance matrix Sigma_epsilon for the measurement error model.
#' @param B Integer. The number of simulated datasets.
#' @param Lambda Numeric vector. A sequence of lambda values for simulated datasets.
#' @param extrapolation Character. A regression model used for extrapolation ("linear", "quadratic", "nonlinear").
#' @param bootstrap_number Numeric. The number of bootstrap samples (default is 250).
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{A data frame with the following columns:
#'     \itemize{
#'       \item \code{Naive_ATE}: Naive estimate of the ATE.
#'       \item \code{Sigma_epsilon}: The covariance matrix Sigma_epsilon for the measurement error model.
#'       \item \code{p10}: The probability of misclassified Y given Y = 0.
#'       \item \code{p11}: The probability of correctly classified Y given Y = 1.
#'       \item \code{Extrapolation}:  A regression model used for extrapolation ("linear", "quadratic", "nonlinear").
#'       \item \code{ATE}: Mean ATE estimate from the bootstrap samples.
#'       \item \code{SE}: Standard error of the ATE estimate.
#'       \item \code{CI}: 95% confidence interval for the ATE estimate.
#'     }
#'   }
#'   \item{boxplot}{A ggplot object representing the boxplot of the ATE estimates.}
#' }
#'
#' @details
#' The `ATE.ERROR.XY` function is designed to handle measurement error in covariates and misclassification in outcomes by using the augmented simulation-extrapolation approach. 
#'
#' @examples
#' \donttest{
#' library(ATE.ERROR)
#' data(Simulated_data)
#' Y_star <- Simulated_data$Y_star
#' A <- Simulated_data$T
#' Z <- Simulated_data$Z
#' X_star <- Simulated_data$X_star
#' p11 <- 0.8
#' p10 <- 0.2
#' sigma_epsilon <- 0.1
#' B <- 100
#' Lambda <- seq(0, 2, by = 0.5)
#' bootstrap_number <- 10
#' result <- ATE.ERROR.XY(Y_star, A, Z, X_star, p11, p10, sigma_epsilon, B, Lambda, 
#'                        "linear", bootstrap_number)
#' print(result$summary)
#' print(result$boxplot)
#' }
#'
#' @import ggplot2
#' @importFrom stats glm predict quantile
#' @export
ATE.ERROR.XY <- function(Y_star, A, Z, X_star, p11, p10, sigma_epsilon, B = 100, 
                         Lambda = seq(0, 2, by = 0.5), extrapolation = "linear", 
                         bootstrap_number = 250) {
  
  Naive_ATE_XY <- Naive_Estimation(Y_star, A, Z, X_star)
  n <- length(Y_star)
  
  # Check if sigma_epsilon is a scalar or a matrix
  if (is.matrix(sigma_epsilon)) {
    # Covariance matrix case (Vector/Matrix inputs)
    Sigma_epsilon <- sigma_epsilon
    X_star_sim_list <- lapply(Lambda, function(lambda) {
      replicate(B, X_star + sqrt(lambda) * mvtnorm::rmvnorm(n, sigma = Sigma_epsilon))
    })
  } else if (is.numeric(sigma_epsilon) && length(sigma_epsilon) == 1) {
    # Scalar case
    X_star_sim_list <- lapply(Lambda, function(lambda) {
      replicate(B, X_star + sqrt(lambda) * rnorm(n, 0, sigma_epsilon))
    })
  } else {
    stop("sigma_epsilon must be either a scalar or a covariance matrix.")
  }
  
  # Step 2: Estimation
  tau_hat_lambda_list <- lapply(1:bootstrap_number, function(iter) {
    tau_hat_lambda <- numeric(length(Lambda))
    for (l in 1:length(Lambda)) {
      X_star_sim <- X_star_sim_list[[l]]
      tau_hat_lambda[l] <- mean(apply(X_star_sim, 2, function(X_sim) {
        sample_idx <- sample(1:n, n, replace = TRUE)
        Y_boot <- Y_star[sample_idx]
        A_boot <- A[sample_idx]
        Z_boot <- Z[sample_idx]
        X_boot <- X_sim[sample_idx]
        logit_model_boot <- glm(A_boot ~ Z_boot + X_boot, family = binomial(link = "logit"))
        e_hat_boot <- predict(logit_model_boot, type = "response")
        E_Y1_boot <- (A_boot * Y_boot / e_hat_boot - p10) / (p11 - p10)
        E_Y0_boot <- ((1 - A_boot) * Y_boot / (1 - e_hat_boot) - p10) / (p11 - p10)
        mean(E_Y1_boot) - mean(E_Y0_boot)
      }))
    }
    tau_hat_lambda
  })
  
  # Step 3: Extrapolation
  tau_tilde_list <- sapply(tau_hat_lambda_list, function(tau_hat_lambda) {
    if (extrapolation == "linear") {
      model <- lm(tau_hat_lambda ~ Lambda)
      tau_tilde <- predict(model, newdata = data.frame(Lambda = -1))
    } else if (extrapolation == "quadratic") {
      model <- lm(tau_hat_lambda ~ poly(Lambda, 2))
      tau_tilde <- predict(model, newdata = data.frame(Lambda = -1))
    } else if (extrapolation == "nonlinear") {
      model <- lm(tau_hat_lambda ~ Lambda / (1 + Lambda))
      tau_tilde <- predict(model, newdata = data.frame(Lambda = -1))
    }
    return(tau_tilde)
  })
  
  # Calculate mean, standard error, and confidence interval
  mean_estimate <- mean(tau_tilde_list)
  se_estimate <- sd(tau_tilde_list)
  ci_estimate <- quantile(tau_tilde_list, probs = c(0.025, 0.975))
  
  # Create the table
  ATE.ERROR.XY_table <- data.frame(
    Naive_ATE_XY = round(Naive_ATE_XY, 3),
    Sigma_epsilon = ifelse(is.matrix(sigma_epsilon), "Covariance Matrix", sigma_epsilon),
    p10 = p10,
    p11 = p11,
    Extrapolation = extrapolation,
    ATE = round(mean_estimate, 3),
    SE = round(se_estimate, 3),
    CI = paste0("(", round(ci_estimate[1], 3), ", ", round(ci_estimate[2], 3), ")")
  )
  
  # Create the boxplot
  ATE.ERROR.XY_data <- data.frame(
    ATE = tau_tilde_list,
    Method = factor(rep(paste("ATE.ERROR.XY", extrapolation, sep = "_"), length(tau_tilde_list)))
  )
  
  list(
    summary = ATE.ERROR.XY_table,
    boxplot = ggplot(ATE.ERROR.XY_data, aes(x = .data$Method, y = .data$ATE, fill = .data$Method)) +
      geom_boxplot() +
      geom_hline(aes(yintercept = Naive_ATE_XY, color = "Naive ATE"), linetype = "dashed") +
      labs(title = paste("ATE Estimates from ATE.ERROR.XY Method (", extrapolation, ")", sep = ""), 
           y = "ATE Estimate") +
      scale_color_manual(name = NULL, values = c("Naive ATE" = "red")) +
      theme_minimal() +
      theme(legend.position = "right") +
      guides(fill = guide_legend(title = NULL, order = 1),
             color = guide_legend(title = NULL, override.aes = list(linetype = "dashed"), order = 2))
    
  )
}
