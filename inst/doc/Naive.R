## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ATE.ERROR)
set.seed(1)
data(Simulated_data)
Y_star <- Simulated_data$Y_star
A <- Simulated_data$T
Z <- Simulated_data$Z
X_star <- Simulated_data$X_star

## -----------------------------------------------------------------------------
Naive_ATE_XY <- Naive_Estimation(Y_star, A, Z, X_star)
print(Naive_ATE_XY)

