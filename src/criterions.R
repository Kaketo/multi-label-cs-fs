source("./src/utils.R")

# SIG with neigborhood
delta_optim <- function(X, gamma){
  sd_vector <- apply(scale(X), 2, sd) / gamma
  return(1/gamma * sum(sd_vector))
}

neighborhood_matrix <- function(X, gamma){
  return(as.matrix(dist(scale(X))) < delta_optim(X, gamma))
}

description_degree <- function(Y, n_mat, index){
  neighborhood_idxs <- n_mat[index,]
  return(colSums(Y[neighborhood_idxs,]) / sum(Y[neighborhood_idxs,]))
}

SIG <- function(X, Y, gamma){
  n_mat <- neighborhood_matrix(X, gamma)
  sig_raw <- t(sapply(1:nrow(X), function(i) description_degree(Y, n_mat, i)))
  return(colMeans(sig_raw))
}

AMI <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  while(length(selected_features) < k){
    if(length(selected_features) == 1){
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2 <- mi_matrix[remaining_features, selected_features]
      J = S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }else if(length(remaining_features) == 1){
      S1 <- sum(mi_matrix[targets, remaining_features])
      S2 <- sum(mi_matrix[remaining_features, selected_features])
      J = S1 - S2
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J)
      )
    }else{
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2 <- rowSums(mi_matrix[remaining_features, selected_features])
      J = S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}

SCLS <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  h_matrix <- hScores(X)
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  while(length(selected_features) < k){
    if(length(selected_features) == 1){
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2_1 <- mi_matrix[remaining_features, selected_features] / h_matrix[remaining_features]
      S2 <- S2_1 * S1
      J = S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }else if(length(remaining_features) == 1){
      S1 <- sum(mi_matrix[targets, remaining_features])
      S2_1 <- sum(mi_matrix[remaining_features, selected_features]) / h_matrix[remaining_features]
      S2 <- S2_1 * S1
      J = S1 - S2
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J)
      )
    }else{
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      if(length(selected_features) == 0){}
      S2_1 <- rowSums(mi_matrix[remaining_features, selected_features]) / h_matrix[remaining_features]
      S2 <- S2_1 * S1
      J = S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}

MDMR <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  while(length(selected_features) < k){
    if(length(selected_features) == 1){
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2_1 <- mi_matrix[remaining_features, selected_features]
      S2_2 <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        for(Yl in targets){
          S2_2 <- S2_2 + cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
        }
      }
      S2 <- S2_1 - S2_2
      J = S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }else if(length(remaining_features) == 1){
      S1 <- sum(mi_matrix[targets, remaining_features])
      S2_1 <- sum(mi_matrix[remaining_features, selected_features])
      S2_2 <- 0
      for(Xj in selected_features){
        for(Yl in targets){
          S2_2 <- S2_2 + cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
        }
      }
      S2 <- S2_1 - S2_2
      J = length(selected_features)*S1 - S2
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J)
      )
    }else{
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2_1 <- rowSums(mi_matrix[remaining_features, selected_features])
      S2_2 <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        for(Yl in targets){
          S2_2 <- S2_2 + cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
        }
      }
      S2 <- S2_1 - S2_2
      ifelse(length(selected_features) == 0 , 1, length(selected_features))
      if(length(selected_features)){}
      J = length(selected_features)*S1 - S2
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}

MVML <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  while(length(selected_features) < k){
    if(length(selected_features) == 1){
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2 <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        for(Yl in targets){
          S2 <- S2 + iiScores(X[,remaining_features], X[,Xj], Y[,Yl])
        }
      }
      S3 <- rep(0,length(remaining_features))
      for(Yl1 in targets){
        for(Yl2 in targets){
          if(Yl1 != Yl2){
            S3 <- S3 + iiScores(X[,remaining_features], Y[,Yl1], Y[,Yl2])
          }
        }
      }
      J = S1 + S2 + S3
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }else if(length(remaining_features) == 1){
      S1 <- sum(mi_matrix[targets, remaining_features])
      S2 <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        for(Yl in targets){
          S2 <- S2 + iiScores(X[,remaining_features], X[,Xj], Y[,Yl])
        }
      }
      S3 <- 0
      for(Yl1 in targets){
        for(Yl2 in targets){
          if(Yl1 != Yl2){
            S3 <- S3 + iiScores(X[,remaining_features], Y[,Yl1], Y[,Yl2])
          }
        }
      }
      J = S1 + S2 + S3
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J)
      )
    }else{
      S1 <- rowSums(mi_matrix[remaining_features, targets])
      S2 <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        for(Yl in targets){
          S2 <- S2 + iiScores(X[,remaining_features], X[,Xj], Y[,Yl])
        }
      }
      S3 <- rep(0,length(remaining_features))
      for(Yl1 in targets){
        for(Yl2 in targets){
          if(Yl1 != Yl2){
            S3 <- S3 + iiScores(X[,remaining_features], Y[,Yl1], Y[,Yl2])
          }
        }
      }
      J = S1 + S2 + S3
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J,
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}

CFSM <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  normalized_costs <- normalize_min_max(costs, 0.1, 0.9)
  normalized_named_costs <- create_named_costs(normalized_costs, X)
  
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  Js <- c()
  
  gamma <- 2
  S1 <- SIG(X, Y, gamma)
  
  while(length(selected_features) < k){
    if(length(remaining_features) == 1){
      S2 <- mi_matrix[remaining_features, targets]
      J = S1 %*% S2
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J * (1-normalized_named_costs[remaining_features]))
      )
    }else{
      S2 <- mi_matrix[remaining_features, targets]
      J = rowSums(S1 * S2)
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J * (1-normalized_named_costs[remaining_features]),
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores, raw_Js=Js))
    }
    
    new_feature <- setdiff(selection_restult$selected_features, selected_features)
    Js <- append(Js, J[new_feature])
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores, raw_Js=Js))
}

PROPOSED_CRITERION <- function(X, Y, costs, budget, lambda){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  Js <- c()
  
  while(length(selected_features) < k){
    if(length(selected_features) == 0){
      mi_matrix <- miMatrix(cbind(X, Y))
      J = (rowSums(mi_matrix[remaining_features, targets]) / length(targets))
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J - lambda * named_costs[remaining_features],
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }else if(length(remaining_features) == 1){
      S <- 0
      for(Xj in selected_features){
        for(Yl in targets){
          S_1 <- cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
          S_2_1 <- rep(chScores(Y[,Yl], X[,Xj]), length(remaining_features))
          S_2_2 <- chScores(X[,remaining_features], X[,Xj])
          S_2 <- rowMins(cbind(S_2_1,S_2_2), value = TRUE)
          
          S <- S + S_1/S_2
        }
      }
      J = (1/(length(selected_features)*length(targets))*S)
      
      selection_restult = list(
        remaining_features = c(),
        selected_features = append(selected_features, remaining_features),
        selected_scores = append(selected_scores, J - lambda * named_costs[remaining_features])
      )
    }else{
      S <- rep(0,length(remaining_features))
      for(Xj in selected_features){
        S1 <- rowSums(sapply(Y, function(Yl) cmiScores(X[,remaining_features],Yl,X[,Xj])))
        M1 <- chScores(Y, X[,Xj])
        M1[M1==0] <- min(M1[M1>0])
        # M2 <- chScores(X[,remaining_features], X[,Xj])
        # M2[M2==0] <- min(M2[M2>0])
        # S2 <- sapply(M2, function(x) min(x,M1))
        S2 <- min(M1)
        
        S <- S + S1/S2
      }
      J = (1/(length(selected_features)*length(targets))*S)
      
      selection_restult <- select_features_based_on_criterion(
        criterion_value = J - lambda * named_costs[remaining_features],
        remaining_features = remaining_features,
        selected_features = selected_features,
        selected_scores = selected_scores
      )
    }
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores, raw_Js=Js))
    }
    
    new_feature <- setdiff(selection_restult$selected_features, selected_features)
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
    
    if(length(selection_restult$remaining_features) == 0){
      Js <- append(Js, J)
    }else{
      Js <- append(Js, J[new_feature])
    }
  }
  
  
  return(list(selection=selected_features, score=selected_scores, raw_Js=Js))
}

STFS <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  q <- ncol(Y)
  
  mi_matrix <- miMatrix(cbind(X, Y))
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  A1 <- rowSums(mi_matrix[colnames(X), targets])
  A1 <- A1 / ncol(Y)
  
  A2 <- rep(0,ncol(X))
  for(Yl in targets){
    for(Yk in targets){
      if(Yl != Yk){
        A2 <- A2 + cmiScores(X, Y[,Yl], Y[,Yk])
      }
    }
  }
  A2 <- A2 / ncol(Y)
  
  while(length(selected_features) < k){
    B1 <- rep(0,length(remaining_features))
    
    for(Xj in selected_features){
      for(Yl in targets){
        B1 <- B1 + cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
      }
    }
    
    B1 <- B1 / (length(selected_features) * ncol(Y))
    
    J <- A1[remaining_features] + A2[remaining_features] + B1
    
    selection_restult <- select_features_based_on_criterion(
      criterion_value = J,
      remaining_features = remaining_features,
      selected_features = selected_features,
      selected_scores = selected_scores
    )
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}

DCRMFS <- function(X, Y, costs, budget){
  X <- data.frame(X)
  Y <- data.frame(Y)
  k <- ncol(X)
  q <- ncol(Y)
  
  mi_matrix <- miMatrix(X)
  named_costs <- create_named_costs(costs, X)
  
  targets <- colnames(Y)
  remaining_features <- colnames(X)
  selected_features <- c()
  selected_scores <- c()
  
  A1 <- rep(0,ncol(X))
  for(Yl in targets){
    for(Yk in targets){
      if(Yl != Yk){
        A1 <- A1 + cmiScores(X, Y[,Yl], Y[,Yk])
      }
    }
  }
  
  while(length(selected_features) < k){
    if(length(remaining_features) == 1){
      A2 <- sum(mi_matrix[remaining_features, selected_features])
    }else if(length(selected_features) == 1){
      A2 <- mi_matrix[remaining_features, selected_features]
    }else{
      A2 <- rowSums(mi_matrix[remaining_features, selected_features])
    }
    B1 <- rep(0,length(remaining_features))
    
    for(Xj in selected_features){
      for(Yl in targets){
        B1 <- B1 + cmiScores(X[,remaining_features], Y[,Yl], X[,Xj])
      }
    }
    
    J <- B1 + A1[remaining_features] - A2
    
    selection_restult <- select_features_based_on_criterion(
      criterion_value = J,
      remaining_features = remaining_features,
      selected_features = selected_features,
      selected_scores = selected_scores
    )
    
    if(sum(named_costs[selection_restult$selected_features]) > budget){
      return(list(selection=selected_features, score=selected_scores))
    }
    
    remaining_features <- selection_restult$remaining_features
    selected_features <- selection_restult$selected_features
    selected_scores <- selection_restult$selected_scores
  }
  
  return(list(selection=selected_features, score=selected_scores))
}
