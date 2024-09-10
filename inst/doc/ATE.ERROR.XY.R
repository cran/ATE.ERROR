## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)

## ----setup--------------------------------------------------------------------
library(ATE.ERROR)
set.seed(1)
data(Simulated_data)
Y_star <- Simulated_data$Y_star
Y <- Simulated_data$Y
A <- Simulated_data$T
Z <- Simulated_data$Z
X_star <- Simulated_data$X_star
X <- Simulated_data$X
p11 <- 0.8
p10 <- 0.2
sigma_epsilon <- 0.1
B <- 100
Lambda <- seq(0, 2, by = 0.5)
bootstrap_number <- 10

## -----------------------------------------------------------------------------
ATE.ERROR.XY_results_linear <- ATE.ERROR.XY(Y_star, A, Z, X_star, p11, p10, sigma_epsilon, 
                                            B, Lambda, extrapolation = "linear", 
                                            bootstrap_number)

## -----------------------------------------------------------------------------
ATE.ERROR.XY_results_quadratic <- ATE.ERROR.XY(Y_star, A, Z, X_star, p11, p10, sigma_epsilon, 
                                               B, Lambda, extrapolation = "quadratic", 
                                               bootstrap_number)


## -----------------------------------------------------------------------------
ATE.ERROR.XY_results_nonlinear <- ATE.ERROR.XY(Y_star, A, Z, X_star, p11, p10, sigma_epsilon, 
                                               B, Lambda, extrapolation = "nonlinear", 
                                               bootstrap_number)

## -----------------------------------------------------------------------------
combined_summary <- rbind(
  ATE.ERROR.XY_results_linear$summary,
  ATE.ERROR.XY_results_quadratic$summary,
  ATE.ERROR.XY_results_nonlinear$summary)

## -----------------------------------------------------------------------------
True_ATE <- True_Estimation(Y, A, Z, X)

## -----------------------------------------------------------------------------
Naive_ATE_XY <- Naive_Estimation(Y_star, A, Z, X_star)

## -----------------------------------------------------------------------------
combined_summary <- data.frame(True_ATE = round(True_ATE, 3), combined_summary)
print(combined_summary)

## ----fig.width=8.5, fig.height=4----------------------------------------------
combined_data <- rbind(
  ATE.ERROR.XY_results_linear$boxplot$data,
  ATE.ERROR.XY_results_quadratic$boxplot$data,
  ATE.ERROR.XY_results_nonlinear$boxplot$data
)

combined_plot <- ggplot(combined_data, aes(x = Method, y = ATE, fill = Method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = Naive_ATE_XY, color = "naive estimate"), 
             linetype = "dashed") +
  geom_hline(aes(yintercept = True_ATE, color = "true estimate"), 
             linetype = "dashed") +
  scale_color_manual(name = NULL, values = c("naive estimate" = "red", 
                                             "true estimate" = "blue")) +
  labs(title = "ATE Estimates from the ATE.ERROR.XY Method with different 
       Approximations of the Extrapolation Function", 
       y = "ATE Estimate") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = NULL, order = 1),
         color = guide_legend(title = NULL, override.aes = list(linetype = "dashed"),
                              order = 2))

print(combined_plot)

