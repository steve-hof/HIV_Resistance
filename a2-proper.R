########## STATS 454 - ASSIGNMENT 2 ##########
### Authors:  Steve Hof and Stephen Scinocca 
###
### Purpose:  Using the HIV data from Assignment 1, implement 8 classification
###           algorithms to classify to judge whether the drug in question is 
###           is resistent to the mutation. To choose the best algorithm we used
###           5 - fold cross validation.
###
### Measurement: Miss-classification rate ( (FP + FN) / 2 )
###
##########
              

rm(list=ls()) # refresh the workspace
setwd("~/school/math/STATS_454/assignments/two/most_recent") # set working directory
set.seed(12) # set random seed so results are reproducible

# Import necessary libraries
library(MASS)
library(caret)
library(glmnet)
library(class)
library(tree)


#########################################
############### FUNCTIONS ###############


########## MODELS ##########
### The models we need to fit are:
###     Binary Y:
###         - Logistic Regression (glm)
###         - Penalized Logistic Regression: LASSO, Elasticnet, Ridge
###         - LDA
###         - KNN
###         - Classification Tree
###
###     Continuous Y:
###         - Linear Regression (from assignment 1)
###         - Penalized Linear Regression: LASSO, Elasticnet, Ridge (from A1)
###         - Regression Tree
##########


########### Model: LDA ##########
###   Input:  Y, X, (train and test sets)
####
###   Output: Matrix of MAP classifications
###########


########### Model: Logistic Regression ##########
###   Input:  Y, xtrain values for this fold, xtest values for this fold
####
###   Output: predicted values on testing data
###########
log_reg <- function(y, x_train, x_test) {
  fit <- glm(y~., data=x_train, family="binomial")
  result <- predict(fit, newdata=x_test, type="response") > 0.5
  return(result)
}


########### Model: LASSO ##########
###   Input:  Y, xtrain values for this fold, xtest values for this fold
####
###   Output: predicted values on testing data
###########
lasso_bin <- function(y, x_train, x_test) {
  fit <- cv.glmnet(data.matrix(x_train[, -1], x_train[, 1], family="binomial"))
  result <- predict(fit1, newx=data.matrix(x_test[, -1]), s="lambda.1se", type="response") > 0.5
  return(result)
}


####################################
############### MAIN ###############

# Import the data
load("HIV.complete.rda")

# Need to reverse the log transform the Y data and we only need to use one drug
# Use 2 as the cutt-off for representing the y values as binary variables (TRUE / FALSE)
y <- 10^YY
y <- y[, 1]
y <- (y > 2)


# Set up matrix as data frame
X <- XX # I prefer X over XX, haha
df <- data.frame(y=y, X)


# Set up folds for cross validation 
num_folds <- 5
fold_idx <- createFolds(y, k=num_folds, list=F)


# Set up matrix to hold all predictions (#obs by #algs)
n <- nrow(df)
algs <- c("Logistic", "Ridge", "LASSO", "Elastic Net", "LDA", "KNN", "Classification Tree",
          "Linear Regression", "Ridge (cont.)", "LASSO (cont.)", "Elastic Net (cont.)",
          "Regression Tree")
predicts <- array(NA, c(n, length(algs)))
dimnames(predicts)[[2]] <- algs

# for(i in 1:num_folds) {
  i = 1
  # set up train / test data (we test on one fold and train on num_folds - 1)
  train_df <- df[fold_idx != i, ]
  test_df <- df[fold_idx == i, ]
  
  # Logistic regression
  predicts[fold_idx == i, c("Logistic")] <- log_reg(y, train_df, test_df)
  
  # Ridge (binary)
  
  # LASSO (binary)
  predicts[fold_idx == i, c("LASSO")] <- lasso_bin(y, train_df, test_df)
  
  # Elastic Net (binary)
  
  # LDA
  
  # KNN
  
  # Classification Tree
  
  # Linear Regression
  
# }