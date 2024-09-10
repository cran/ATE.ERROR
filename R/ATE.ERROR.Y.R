#' ATE.ERROR.Y Function for Estimating Average Treatment Effect (ATE) with Misclassification in Y
#'
#' This function performs estimation of the Average Treatment Effect (ATE) using the ATE.ERROR.Y method,
#' which accounts for misclassification in the binary outcome variable Y. The method calculates consistent estimates
#' of the ATE in the presence of misclassified outcomes by leveraging logistic regression and bootstrap sampling.
#'
#' @param Y_star Numeric vector. The observed binary outcome variable, which may be subject to misclassification.
#' @param A Numeric vector. The binary treatment indicator (1 if treated, 0 if control).
#' @param Z Numeric vector. A precisely measured covariate vector.
#' @param X Numeric vector. A precisely measured covariate vector.
#' @param p11 Numeric. The probability of correctly classified Y given Y = 1.
#' @param p10 Numeric. The probability of misclassified Y given Y = 0.
#' @param bootstrap_number Integer. The number of bootstrap samples (default is 250) used to obtain the associated variance estimate.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{A data frame with the following columns:
#'     \itemize{
#'       \item \strong{Naive_ATE}: Naive estimate of the ATE, ignoring misclassification.
#'       \item \strong{ATE}: Mean ATE estimate from the bootstrap samples, accounting for misclassification.
#'       \item \strong{SE}: Standard error of the ATE estimate.
#'       \item \strong{CI}: 95% confidence interval for the ATE estimate.
#'     }
#'   }
#'   \item{boxplot}{A ggplot object representing the boxplot of the ATE estimates.}
#' }
#'
#' @details
#' The function first calculates consistent estimates of the ATE, correcting for misclassification in the outcome variable Y.
#' The logistic model is used to estimate the propensity scores for the treatment assignment, which are then adjusted using the provided
#' misclassification probabilities p11 and p10. Bootstrap sampling is performed to estimate the variance and construct confidence intervals
#' for the ATE estimates.
#'
#' @examples
#' library(ATE.ERROR)
#' data(Simulated_data)
#' Y_star <- Simulated_data$Y_star
#' A <- Simulated_data$T
#' Z <- Simulated_data$Z
#' X <- Simulated_data$X
#' p11 <- 0.8
#' p10 <- 0.2
#' bootstrap_number <- 250
#' result <- ATE.ERROR.Y(Y_star, A, Z, X, p11, p10, bootstrap_number)
#' print(result$summary)
#' print(result$boxplot)
#'
#' @export
ATE.ERROR.Y <- function(Y_star, A, Z, X, p11, p10, bootstrap_number = 250) {
  Naive_ATE_Y <- Naive_Estimation(Y_star, A, Z, X)
  n <- length(Y_star)
  estimates <- replicate(bootstrap_number, {
    sample_idx <- sample(1:n, n, replace = TRUE)
    Y_boot <- Y_star[sample_idx]
    A_boot <- A[sample_idx]
    Z_boot <- Z[sample_idx]
    X_boot <- X[sample_idx]
    
    # Logistic model to estimate propensity scores
    logit_model_boot <- glm(A_boot ~ Z_boot + X_boot, family = binomial(link = "logit"))
    e_hat_boot <- predict(logit_model_boot, type = "response")
    
    # Estimating E^*(Y(1)) and E^*(Y(0))
    E_Y1_boot <- (A_boot * Y_boot / e_hat_boot - p10) / (p11 - p10)
    E_Y0_boot <- ((1 - A_boot) * Y_boot / (1 - e_hat_boot) - p10) / (p11 - p10)
    
    # Calculating the ATE
    mean(E_Y1_boot) - mean(E_Y0_boot)
  })
  
  # Calculate mean, standard error, and confidence interval
  mean_estimate <- mean(estimates)
  se_estimate <- sd(estimates)
  ci_estimate <- quantile(estimates, probs = c(0.025, 0.975))
  
  # Create the table
  ATE.ERROR.Y_table <- data.frame(
    Naive_ATE_Y = round(Naive_ATE_Y, 3),
    ATE = round(mean_estimate, 3),
    SE = round(se_estimate, 3),
    CI = paste0("(", round(ci_estimate[1], 3), ", ", round(ci_estimate[2], 3), ")")
  )
  
  # Create the boxplot
  ATE.ERROR.Y_data <- data.frame(
    ATE = estimates,
    Method = factor(rep("ATE.ERROR.Y", length(estimates)))
  )
  
  list(
    summary = ATE.ERROR.Y_table,
    boxplot =  ggplot(ATE.ERROR.Y_data, aes(x = .data$Method, y = .data$ATE, fill = .data$Method)) +
      geom_boxplot() +
      geom_hline(aes(yintercept = Naive_ATE_Y, color = "naive estimate"), linetype = "dashed") +
      labs(title = "ATE Estimates from ATE.ERROR.Y Method", y = "ATE Estimate") +
      scale_color_manual(name = NULL, values = c("naive estimate" = "red")) +
      theme_minimal() +
      theme(legend.position = "right") +
      guides(fill = guide_legend(title = NULL, order = 1),
             color = guide_legend(title = NULL, override.aes = list(linetype = "dashed"), order = 2))
  )
}
