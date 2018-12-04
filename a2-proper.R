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
setwd("~/school/math/STATS_454/assignments/two/working_files") # set working directory
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
###


########### Model: LDA ##########
###   Input:  Y, X, (train and test sets)
####
###   Output: Matrix of MAP classifications
###########
lda_func <- function(y, x_train, x_test) {
  fit <- lda(y~., data=x_train)
  
  result <- predict(fit, newdata=x_test)$class #- 1
  return(result)
}

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
lasso_bin <- function(y, train_temp, test_temp, alpha) {
  x_train <- data.matrix(train_temp)
  x_test <- data.matrix(test_temp)
  fit <- cv.glmnet(x_train[, -1], x_train[, 1], family="binomial", alpha=alpha)
  result <- predict(fit, newx=x_test[, -1], s="lambda.1se", type="response") > 0.5
  return(result)
}


########### Model: KNN ##########
###   Input:  Y, xtrain values for this fold, xtest values for this fold
####
###   Output: predicted values on testing data
####
###   KNN doesn't have a built in option for cross-validation 
###   (to choose optimum k value), so this function adds in that
###   necessary cross validation loop
###########
knn_func <- function(y, x_train, x_test) {
  num_rows <- nrow(x_train)
  num_trials <- 4
  num_knn_folds <- 10
  
  # Creat index of folds for comparing k's in cross validation
  knn_fold_idx <- createFolds(x_train$y, k=num_knn_folds, list=F)
  
  # Set up empty matrix to hold results for each round of testing number of folds to use
  k_results <- matrix(NA, num_rows, num_trials)
  
  # Perform cross validation to maximize training potential
  for(curr_fold in 1:num_knn_folds) {
    k_train <- x_train[knn_fold_idx != curr_fold, ]
    k_test <- x_train[knn_fold_idx == curr_fold, ]
    
    # Loop through and fill out k_results matrix
    for(curr_k in 1:num_trials) {
      k_results[knn_fold_idx == curr_fold, curr_k] <- (k_test[, 1] != knn(k_train[, -1], k_test[, -1], 
                                                                          k_train[, 1], k=curr_k))
    }
  }
  
  # Count up number of errors in each column
  totals <- colSums(k_results)
  
  # Choose column with fewest errors (tie goes to smaller k)
  k <- 1
  low <- totals[1]
  for (i in 1:length(totals)) {
    if (totals[i] < low) {
      low <- totals[i]
      k <- i
    }
  }
  
  result <- knn(x_train[, -1], x_test[, -1], x_train[, 1], k=k)
  return(result)
}


########### Model: Classification Tree ##########
####
###   Purpose: Use Tree library to train, prune, and test HIV data
####
###   Input:  Y, xtrain values for this fold, xtest values for this fold
####
###   Output: predicted values on testing data
###########
class_tree <- function(y, x_train, x_test) {
  dummy_tree <- tree(factor(y)~., x_train)
  dummy_cv <- cv.tree(dummy_tree, FUN=prune.misclass)
  tree_size <- dummy_cv$size[which.min(dummy_cv$dev)]
  fit <- prune.misclass(dummy_tree, best=tree_size)
  result <- predict(fit, x_test, type="class")
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


# Set up matrix to hold all predictions (dimension: #obs X #algs)
n <- nrow(df)
algs <- c("Logistic", "Ridge", "LASSO", "Elastic Net", "LDA", "KNN", "Classification Tree",
          "Linear Regression", "Ridge (cont.)", "LASSO (cont.)", "Elastic Net (cont.)",
          "Regression Tree")
predicts <- array(NA, c(n, length(algs)))
dimnames(predicts)[[2]] <- algs

# for(i in 1:num_folds) {
  curr_fold = 1
  # set up train / test data (we test on one fold and train on num_folds - 1)
  train_df <- df[fold_idx != i, ]
  test_df <- df[fold_idx == i, ]
  
  # Logistic regression
  cat("fitting Logistic Regression (glm) for fold", curr_fold, "...\n")
  # predicts[fold_idx == i, c("Logistic")] <- log_reg(y, train_df, test_df)
  
  # Ridge (binary)
  cat("fitting Ridge (binary) for fold", curr_fold, "...\n")
  # predicts[fold_idx == i, c("Ridge")] <- lasso_bin(y, train_df, test_df, alpha=0)
  
  # Elastic Net
  cat("fitting Elastic Net (binary) for fold", curr_fold, "...\n")
  # predicts[fold_idx == i, c("Elastic Net")] <- lasso_bin(y, train_df, test_df, alpha=.5)
  
  # LASSO (binary)
  cat("fitting LASSO (binary) for fold", curr_fold, "...\n")
  # predicts[fold_idx == i, c("LASSO")] <- lasso_bin(y, train_df, test_df, alpha=1)
  
  # LDA
  cat("fitting LDA for fold", curr_fold, "...\n")
  predicts[fold_idx == curr_fold, c("LDA")] <- lda_func(y, train_df, test_df)
  
  # KNN
  cat("Fitting KNN for fold", curr_fold, "...\n")
  predicts[fold_idx == curr_fold, c("KNN")] <- knn_func(y, train_df, test_df)
  
  # Classification Tree
  cat("Fitting Classification Tree for fold", curr_fold, "...\n")
  predicts[fold_idx == curr_fold, c("Classification Tree")] <- class_tree(y, train_df, test_df)
  
  
  

  # Linear Regression
  
# }