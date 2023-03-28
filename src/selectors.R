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