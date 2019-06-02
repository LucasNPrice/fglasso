# Luke Price
# 2019.05.16
# Determines edges from fglasso precision matrix 


getEdges <- function(blocked_matrix) {
  
  if (dim(blocked_matrix)[1] != dim(blocked_matrix)[2]) {
    warning("Blocked matrix is not square!")
  }
  
  p <- dim(blocked_matrix)[1]
  edges <- matrix(0, nrow = p, ncol = p)
  
  for (row in 1:p) {
    for (col in 1:p) {
      
      if (all(blocked_matrix[[row,col]] == 0)) {
        edges[row,col] <- 0
      } else {
        edges[row,col] <- 1
      }
    }
  }
  return(edges)
}

#-------------------------------------------------

getErrors <- function(estimated_edges, simulated_edges) {
  
  conf_mat <- table(estimated_edges, simulated_edges)
  TP <- conf_mat[2,2]
  FP <- conf_mat[2,1]
  TN <- conf_mat[1,1]
  FN <- conf_mat[1,2]
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  return(list("TPR" = TPR, "FPR" = FPR))
}




