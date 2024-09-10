## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)

## -----------------------------------------------------------------------------
library(ATE.ERROR)
set.seed(1)
data(Simulated_data)
Y_star <- Simulated_data$Y_star
A <- Simulated_data$T
Z <- Simulated_data$Z
X_star <- Simulated_data$X_star
X <- Simulated_data$X
Y <- Simulated_data$Y
p11 <- 0.8
p10 <- 0.2
bootstrap_number <- 1000

## -----------------------------------------------------------------------------
set.seed(1)
result <- ATE.ERROR.Y(Y_star, A, Z, X, p11, p10, bootstrap_number)

## -----------------------------------------------------------------------------
True_ATE <- True_Estimation(Y, A, Z, X)

result_summary <- result$summary
result_summary <- data.frame(True_ATE = round(True_ATE, 3), result_summary)
print(result_summary)

## ----fig.width=6, fig.height=4, comment = NA----------------------------------
boxplot_with_true <- result$boxplot +
  geom_hline(aes(yintercept = True_ATE, color = "true estimate"), 
             linetype = "dashed") +
  scale_color_manual(name = NULL, values = c("naive estimate" = "red", 
                                             "true estimate" = "blue")) +
  labs(title = "ATE Estimates from ATE.ERROR.Y Method", y = "ATE Estimate") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = NULL, order = 1),
         color = guide_legend(title = NULL, override.aes = list(linetype = "dashed")
                              , order = 2))

print(boxplot_with_true)

