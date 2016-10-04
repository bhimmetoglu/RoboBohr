# Burak Himmetoglu
# 09-29-2016
#
#
# Analyze GDB data from NJP paper: NJP 15 (2013) 095003 
#
# Libraries
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/home/burak/Works/RoboBohr/GDB")

# Read outcomes 
colnamesOut <- c("AE", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "HOMO", "LUMO", "gwHOMO", "gwLUMO", "Alpha", "scsAlpha")
outcomes <- read.csv("outcomes.csv", header = FALSE )
colnames(outcomes) <- colnamesOut
Y <- outcomes$AE # We are interested in PBE0 atomization energies only

# Convert to eV
eV_to_kcalMol <- 23.0609
Y <- Y/eV_to_kcalMol

# Read Columb Matrix eigenvalues
nDim <- 23
colnamesCM <- paste0("px", 1:nDim)
matrices <- read.csv("coulombL.csv", header = FALSE)
colnames(matrices) <- colnamesCM

# Scale factor
scl <- -max(abs(Y))

# Mean variation
sqrt(mean((mean(Y)-Y)^2))*eV_to_kcalMol

# Split into training and testing
set.seed(101)
inTrain <- sample(1:dim(matrices)[1], size = floor(0.7*dim(matrices)[1]), replace = FALSE)
train.Y <- Y[inTrain]; test.Y <- Y[-inTrain]
train.X <- matrices[inTrain,]; test.X <- matrices[-inTrain,]

## Check linear model:
# Train
library(glmnet)
mod.glm <- cv.glmnet(x = as.matrix(train.X), y = train.Y, nfolds = 5, type.measure = "mse")
pred <- predict(mod.glm, newx = as.matrix(test.X))

# 5-fold CV error and test error
min(sqrt(mod.glm$cvm))*eV_to_kcalMol # This is 5-fold CV error
sqrt(mean((pred - test.Y)^2))*eV_to_kcalMol # This is test error

## XGBoost
# XGBoost style matrices
library(xgboost)
dtrain.X <- xgb.DMatrix(as.matrix(train.X), label = train.Y/scl)
dtest.X <- xgb.DMatrix(as.matrix(test.X), label = test.Y/scl)

# Watchlist
watchlist <- list(train=dtrain.X, test=dtest.X)

# Train by cross-validation
library(caret)
xgb_grid <- expand.grid(eta = 2^seq(-6,-4), colsample_bytree = c(0.2,0.4,0.6),
                        max_depth = c(2,6,8,16), min_child_weight = c(2,6,8,10), 
                        gamma = c(0,1e-4,0.001,0.01))

# Loop over parameters
cv.nround = 1000
cv.nfold = 5
cv.results <- data.frame(nrounds = NULL, rmse = NULL)
for (ind in 1:dim(xgb_grid)[1]){
  # Parameters
  param <- list(booster="gbtree",
                eval_metric="rmse",
                eta = xgb_grid[ind, 1],
                colsample_bytree = xgb_grid[ind,2],
                max_depth = xgb_grid[ind, 3],
                min_child_weight = xgb_grid[ind, 4],
                gamma = xgb_grid[ind, 5],
                lambda = xgb_grid[ind, 6],
                subsample = 0.8)
  # Cross validation
  mod.cv <- xgb.cv(data = dtrain.X, params = param, nfold = cv.nfold, nrounds = cv.nround,
                   early.stop.round = 3, verbose = FALSE)
  
  # Save results
  nround_final <- dim(mod.cv)[1]
  cv.results[ind,1] <- nround_final; cv.results[ind,2] <- mod.cv[nround_final]$test.rmse.mean
  cat("Trained ", ind, " of ", dim(xgb_grid)[1], " :  rmse= ", cv.results[ind,2],  "\n")
}

save(cv.results, file = "cvResultsL.RData")

ind.best <- which.min(cv.results$V2)
nrounds_fit <- cv.results[ind.best, 1] 
# Best parameters (Reference):
# nrounds_fit = 797
#        eta colsample_bytree max_depth min_child_weight gamma
#  0.03125              0.6         6                6     0

###### Final model

# Parameters (trained above) 
param <- list(booster="gbtree",
              eval_metric="rmse",
              eta=0.03125,
              colsample_bytree = 0.6,
              max_depth = 6,
              min_child_weight = 6,
              gamma = 0.0,
              lambda = 1.0,
              subsample = 0.8)

# Test xgboost
xgb.model <- xgb.train(data=dtrain.X, params = param, watchlist=watchlist,nround = nrounds_fit)

# Predict
pred <- predict(xgb.model, newdata = dtest.X)*scl

# RMSE
sqrt(mean((pred - test.Y)^2))*eV_to_kcalMol # 14.43097 kcal/mol
