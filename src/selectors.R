source("./src/criterions.R")
source("./src/utils.R")

lambda_max <- function(X, Y, costs, budget){
  # TODO: Why there is NaN in score at the end of vector in BIBTEX dataset
  l0r <- PROPOSED_1(X, Y, costs, budget=budget, lambda=0)
  J_max <- max(l0r$score)
  
  sc <- sort(unique(costs))
  l_sc <- length(sc)
  mcd <- abs(sc[2] - sc[1])
  if(l_sc == 1){
    stop("All costs are equal to the same value.")
  }else if(l_sc == 2){
    return(J_max / mcd)
  }else{
    for(i in 2:(l_sc - 1)){
      a <- sc[i]
      b <- sc[i+1]
      if(abs(a-b) > 0){
        mcd <- min(mcd, abs(a-b))
      }
    }
  }
  return(J_max / mcd)
}

lambda_optimal <- function(X, Y, costs, budget, l_cnt=10, l_max=NULL){
  l_min <- 0
  if(is.null(l_max)){
    l_max <- lambda_max(X, Y, costs, budget)
  }
  
  l_seq <- seq(l_min, l_max, length.out = l_cnt + 1)[-1]
  l_Js_sum <- c()
  for(l in l_seq){
    J_results <- PROPOSED_1(X, Y, costs, budget, l)
    l_Js_sum <- append(l_Js_sum, sum(J_results$raw_Js))
  }
  argmax <- which.max(l_Js_sum)
  return(l_seq[argmax])  
}


select_features_based_on_criterion <- function(criterion_value, remaining_features, selected_features, selected_scores){
  argmax <- which.max(ifelse(is.na(criterion_value), 0, criterion_value))
  argmax_feature <- names(argmax)
  argmax_score <- max(ifelse(is.na(criterion_value), 0, criterion_value))
  
  
  remaining_features <- setdiff(remaining_features, argmax_feature)
  selected_features <- append(selected_features, argmax_feature)
  selected_scores <- append(selected_scores, argmax_score)
  
  return(list(
    remaining_features = remaining_features,
    selected_features = selected_features,
    selected_scores = selected_scores
  ))
}

score_final_feature_selection_revision <- function(mldr_dataset_train, mldr_dataset_test, costs, budget, lambda){
  X_train <- mldr_dataset_train$dataset[mldr_dataset_train$attributesIndexes]
  Y_train <- mldr_dataset_train$dataset[-mldr_dataset_train$attributesIndexes]
  
  named_costs <- create_named_costs(costs, X_train)
  start <- Sys.time ()
  selection_STFS <- STFS(X_train, Y_train, costs, budget)
  selection_DCRMFS <- DCRMFS(X_train, Y_train, costs, budget)
  selection_stop <- Sys.time ()
  print(paste(" - selection time: ", round(difftime(selection_stop, start, units="secs"), 5), "seconds"))
  
  selected_features <- list()
  selected_features[['STFS']] <- selection_STFS$selection
  selected_features[['DCRMFS']] <- selection_DCRMFS$selection
  results_matrix <- data.frame(matrix(ncol = 27, nrow = 0))
  
  start <- Sys.time ()
  for (name in c("STFS", "DCRMFS")) {
    selected_features_tmp <- selected_features[[name]]
    if(is.null(selected_features_tmp)){
      results_matrix_tmp <- data.frame(
        "criterion" = name,
        "lambda" = NA,
        "cost" = 0,
        "budget" = budget,
        "accuracy"= NA,
        "average.precision"=NA,
        "clp"=NA,
        "coverage"=NA, 
        "F1"=NA,
        "hamming.loss"=NA,
        "macro.AUC"=NA,
        "macro.F1"=NA,
        "macro.precision"=NA ,
        "macro.recall"=NA,
        "margin.loss"=NA,
        "micro.AUC"=NA, 
        "micro.F1"=NA,
        "micro.precision"=NA,
        "micro.recall"=NA,
        "mlp"=NA, 
        "one.error"=NA,
        "precision"=NA,
        "ranking.loss"=NA ,
        "recall"=NA, 
        "subset.accuracy"=NA,
        "wlp"=NA,
        "selection_type"="traditional"
      )
    }else{
      mldr_dataset_train_subset <- create_subset(mldr_dataset_train, cols = selected_features_tmp)
      mldr_dataset_test_subset <- create_subset(mldr_dataset_test, cols = selected_features_tmp)
      
      # model <- mlknn(mldr_dataset_train_subset, k=10)
      # preds <- predict(model, mldr_dataset_test_subset)
      model <- br(mldr_dataset_train_subset, 'KNN', k=10)
      preds <- predict(model, mldr_dataset_test_subset)
      
      results_matrix_tmp <- multilabel_evaluate(mldr_dataset_test_subset, preds)
      results_matrix_tmp['cost'] <- sum(named_costs[selected_features_tmp])
      results_matrix_tmp['budget'] <- budget
      results_matrix_tmp["criterion"] <- name
      results_matrix_tmp['lambda'] <- NA
      results_matrix_tmp['selection_type'] <- "traditional"
      
      results_matrix_tmp <- t(as.data.frame(results_matrix_tmp))
      rownames(results_matrix_tmp) <- NULL
      results_matrix_tmp <- as.data.frame(results_matrix_tmp)
      colnames(results_matrix_tmp) <- gsub("-",".",colnames(results_matrix_tmp))
      results_matrix_tmp <- results_matrix_tmp[,c(
        "criterion",
        "lambda",
        "cost",
        "budget",
        "accuracy", "average.precision", "clp",
        "coverage", "F1", "hamming.loss", "macro.AUC","macro.F1", "macro.precision", "macro.recall",
        "margin.loss", "micro.AUC", "micro.F1", "micro.precision","micro.recall", "mlp", "one.error",
        "precision", "ranking.loss", "recall", "subset.accuracy", "wlp", "selection_type")
      ]
    }
    results_matrix <- rbind(results_matrix, results_matrix_tmp)
  }
  selection_stop <- Sys.time ()
  print(paste(" - scoring time:   ", round(difftime(selection_stop, start, units="secs"), 5), "seconds"))
  
  return(results_matrix)
}
