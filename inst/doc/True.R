## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ATE.ERROR)
set.seed(1)
data(Simulated_data)
Y <- Simulated_data$Y
A <- Simulated_data$T
Z <- Simulated_data$Z
X <- Simulated_data$X

## -----------------------------------------------------------------------------
True_ATE <- True_Estimation(Y, A, Z, X)
print(True_ATE)

