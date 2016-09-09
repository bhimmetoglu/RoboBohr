# Burak Himmetoglu
# begin: 09-08-2016
#
#
# Read pixels combine data and explore
library(dplyr)
library(ggplot2)
library(data.table)

# Working dir
setwd("/home/burak/Works/RoboBohr/methodsVersion/")

## Read the data
# Outcomes first
outcomes <- fread("./data/out.dat.0_16274", skip = 2, header = FALSE, colClasses = c("integer", "numeric"))
colnames(outcomes) <- c("Id", "E")

# AESub (piece to subtract from outcomes to get Eat)
AeSub <- fread("./data/AEsub.out", header = FALSE, colClasses = c("integer","numeric"))
colnames(AeSub) <- c("Id", "Esub")

# Merge AeSub and outcomes; compute atomization energies
outcomesAe <- merge(outcomes, AeSub, by = "Id")
outcomesAe[,Eat := E-Esub]
outcomesAe[,E:=NULL]; outcomesAe[,Esub:=NULL] # Remove unnescessary columns
rm(outcomes,AeSub); gc()

# Scale 
scl <- -max(abs(outcomesAe$Eat), na.rm = TRUE)

# Read the Coulomb Matrix (eigenvalues)
CoulombLambda <- fread("./data/coulombL.csv", header = FALSE)

# Fix column names
nam <- paste0('px', 1:(ncol(CoulombLambda)-1))
nam <- c("Id", nam)
colnames(CoulombLambda) <- nam
CoulombLambda[,Id:=as.integer(Id)] # Make Id variable integer

# Match with Id's so that there is no mistmatch in order
combined <- merge(CoulombLambda, outcomesAe, by="Id")
rm(outcomesAe,CoulombLambda); gc()

# Remove NA's
l.complete <- complete.cases(combined$Eat)
combined <- combined[l.complete,]

# Split into training and testing
set.seed(101)
inTrain <- sample(1:dim(combined)[1], size = floor(0.7*dim(combined)[1]), replace = FALSE)
train.Y <- combined[inTrain, ]$Eat; test.Y <- combined[-inTrain, ]$Eat

# Feature matrices
combined[,Eat:=NULL] # No need for Eat 
combined[,Id:=NULL] # No need for Id
train.X <- combined[inTrain,]; test.X <- combined[-inTrain,]
rm(combined,l.complete); gc()

# XGBoost style matrices
library(xgboost)
dtrain.X <- xgb.DMatrix(as.matrix(train.X), label = train.Y/scl)
dtest.X <- xgb.DMatrix(as.matrix(test.X), label = test.Y/scl)

# Watchlist
watchlist <- list(train=dtrain.X, test=dtest.X)

# Parameters
param <- list(booster="gbtree",
              eval_metric="rmse",
              eta=0.0156,
              colsample_bytree = 0.4,
              max_depth = 16,
              min_child_weight = 10,
              gamma = 0.0,
              lambda = 1.0,
              subsample = 0.8)

# Test xgboost
xgb.model <- xgb.train(data=dtrain.X, params = param, watchlist=watchlist,nround = 600)

# Predict
pred <- predict(xgb.model, newdata = dtest.X)*scl

# RMSE
sqrt(mean((pred - test.Y)^2)) # 0.1309823
