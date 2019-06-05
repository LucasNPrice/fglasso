# block_matrix.R
# blocks and unblocks (square) matrices 

block_matrix <- function(inputmatrix, outersize, innersize){
  # outersize = length of outer matrix ( number of functions/nodes represented in full matrix)
  # innersize = length of inner  matrix (number of harmonics representing each function/node)
  inner_block <- list()
  blocked_matrix <- matrix(inner_block, nrow = outersize, ncol = outersize)
  
  for (i in 1:outersize) {
    for (j in 1:outersize) {
      blocked_matrix[[i,j]] <- inputmatrix[c(((i*innersize)-(innersize-1)):(i*innersize)), c(((j*innersize)-(innersize-1)):(j*innersize))]
    }
  }
  return(blocked_matrix)
}

unblock_matrix <- function(blocked_matrix) {
  p <- dim(blocked_matrix)[1]
  tmp_matrix <- vector("list", p)
  
  for (row in 1:p) {
    tmp_matrix[[row]] <- do.call(cbind, blocked_matrix[row,])
  }
  tmp_matrix <- do.call(rbind, tmp_matrix)
  return(tmp_matrix)
}

