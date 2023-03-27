library(praznik)
library(stringr)

normalize_min_max <- function(x, a, b){
  (b - a)*(x - min(x))/(max(x)-min(x))+a
}

iiScores <- function(X,Y,Z){
  # X - matrix
  # Y - vector
  # Z - vector
  # Calculates II for all combinations of X, Y and Z
  ii <- cmiScores(X, Y, Z) - miScores(X,Y)
  return(ii)
}

chScores <- function(X,Y){
  jh <- jhScores(X,Y) - hScores(X)
  return(jh)
}

create_named_costs <- function(costs, X){
  stopifnot(ncol(X) == length(costs))
  named_costs <- structure(costs, names=colnames(X))
  return(named_costs)
}

train_test_split <- function(X,Y,test_size=0.1, seed=123){
  set.seed(seed)
  smp_size <- floor((1-test_size) * nrow(X))
  train_ind <- sample(seq_len(nrow(X)), size = smp_size)
  X_train <- X[train_ind, ]
  X_test <- X[-train_ind, ]
  Y_train <- Y[train_ind,]
  Y_test <- Y[-train_ind,]
  
  rownames(X_train) <- NULL
  rownames(Y_train) <- NULL
  rownames(X_test) <- NULL
  rownames(Y_test) <- NULL
  
  return(list(X_train=X_train, X_test=X_test, Y_train=Y_train, Y_test=Y_test))
}