## STATS 454 - Assignment 2 

rm(list=ls())
setwd("~/school/math/STATS_454/assignments/two")
library(MASS)
library(tree)

#library(ISLR)
#library(scatterplot3d)
#library(robustbase)

## Fitting Models

# create fake data
r <- 100
c <- 10
hundredths <- seq(from=0, to=1, by=.01)

X <- matrix(rbinom(r*c, 1, 0.5), r, c)
colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")

y.bin <- matrix(rbinom(r, 1, 0.6), r, 1)
data.bin <- cbind(y.bin, X)

y.cont <- sample(hundredths, size=25, replace=TRUE)

data.cont <- cbind(y.cont, X)

# data needs to be a data frame
data.cont <- as.data.frame(data.cont)
data.bin <- as.data.frame(data.bin)

# train test split
train.cont <- sample(1:nrow(data.cont), nrow(data.cont)/2)
train.bin <- sample(1:nrow(data.bin), nrow(data.bin)/2)

# Trees (Regression)
set.seed(12)

tree.cont <- tree(y.cont~., data.cont, subset=train.cont)

summary(tree.cont)
plot(tree.cont)
text(tree.cont,pretty=0)

# make predictions
yhat <- predict(tree.cont, newdata=data.cont[-train.cont,])
data.cont.test <- data.cont[-train.cont, "y.cont"]
plot(yhat, data.cont.test)
abline(0,1)
mean((yhat-data.cont.test)^2)

# LDA
lda <- lda(y~., data, subset=train)
lda
