# Luke Price
# 2019.04.07
# Perform fglasso

require(doParallel)
require(pracma)

numCores <- detectCores()
registerDoParallel(cores = numCores - 1)
getDoParWorkers()


############################################################################
###         ___Main___ Algorithm                                         ###
############################################################################
Algorithm_1 <- function(bmatrix, gamma) {
  
  start_time <- Sys.time()
  print(paste("Number Cores in Use ---", getDoParWorkers()))
  
  p <- dim(bmatrix)[1]
  M <- unique(sapply(bmatrix, nrow))
  
  inner_block <- list(matrix(0, nrow = M, ncol = M))
  Sigma <- matrix(inner_block, nrow = p, ncol = p)
  Theta <- matrix(inner_block, nrow = p, ncol = p)
  diag(Sigma) <- list(diag(M))
  diag(Theta) <- list(diag(M))
  Theta_inv <- Theta
  
  #---------------------------------------------------------------#
  
  error_criterion = 0.001
  norm_error <- c(1)
  iteration = 1
  
  while (tail(norm_error, n = 1) > error_criterion) {
    print(paste("Iteration Time -----", Sys.time() - start_time))
    last_sigma <- Sigma
    
    for (j in 1:p) {
      j_time <- Sys.time()
      
      print(paste("Function: ", j, " / ", p))
      
      # print("Entering updateTheta")
      Theta_inv <- updateTheta(Theta = Theta, Sigma = Sigma, Theta_inv = Theta_inv, p = p, j = j)

      # print("Entering Algorithm 3")
      Theta <- Algorithm_3(Theta = Theta, Theta_inv = Theta_inv, SampleCov = bmatrix, gamma = gamma, p = p, M = M, j = j)
      
      # print("Entering get_U")
      Uj <- get_U(Theta_inv = Theta_inv, Theta = Theta, p = p, M = M, j = j)
      
      # print("Entering updateSigma")
      Sigma <- updateSigma(Sigma = Sigma, Theta_inv = Theta_inv, sampleCov = bmatrix, Uj = Uj, p = p, j = j)
      
      print(paste("J time -----", Sys.time() - j_time))
    }
    
    fNorm <- get_Norm(Sigma = Sigma, last_sigma = last_sigma, p = p)
    norm_error[iteration] <- fNorm
    
    print("--------------------------------------------------")
    print("--------------------------------------------------")
    print(paste("ITERATION ---", iteration, "--- ERROR ---", norm_error[iteration]))
    print("--------------------------------------------------")
    print("--------------------------------------------------")
    iteration <- iteration + 1
    
    }
  
  print(paste("Time to Convergence ---", Sys.time() - start_time))
  
  PrecisionMatrix <- get_Precision(Theta = Theta, SampleCov = SampleCov, Theta_inv = Theta_inv, p = p, M = M)
  
  return(list("PrecisionMatrix" = PrecisionMatrix, "sigma" = Sigma, "theta_inv" = Theta_inv, 
              "numIterations" = iteration, "finalNormError" = fNorm, "normErrors" = norm_error))
  
}


############################################################################
###           Step 2(a) of algorithm 1                                   ###
############################################################################
updateTheta <- function(Theta, Sigma, Theta_inv, p, j) {
  
  SigSmallTranspose <- t(matrix(lapply(Sigma[-j,j], t), ncol = length(Sigma[-j,j])))
  
  
  temp_theta <- foreach(row_block = 1:(p-1)) %dopar% {

    foreach(col_block = 1:(p-1)) %do% {

      temp_theta  <- (Sigma[-j,-j][[row_block, col_block]]
                      - Sigma[-j,j][[row_block]] %*% solve(Sigma[[j,j]]) %*% SigSmallTranspose[[col_block]])

    }
  }
  
  
  for (row_block in 1:length(temp_theta)) {
    
    for (col_block in 1:length(temp_theta[[row_block]])) {
      
      Theta_inv[-j,-j][[row_block, col_block]] <- temp_theta[[row_block]][[col_block]]
      
    }
  }

  return(Theta_inv)
}

##################################################################################
###    Algorithm 3: solve for w_j                                              ###
##################################################################################
Algorithm_3 <- function(Theta, Theta_inv, SampleCov, gamma, p, M, j) {
  
  error_criterion = 0.001
  iteration = 1
  error <- c(1)
  
  while (tail(error, n = 1) > error_criterion) {
    
    iter_time <- Sys.time()
    # print(paste("iteration ---", iteration))
    last_w <- Theta[-j,j]
    theta_nj_transpose <- t(matrix(lapply(Theta_inv[-j,-j], t), nrow = nrow(Theta_inv[-j,-j]), ncol = ncol(Theta_inv[-j,-j])))

    w_j <- foreach(k = 1:(p-1)) %dopar% {

      block_resid = matrix(0, nrow = M, ncol = M)

      foreach(l = 1:(p-1)) %do% {

        if (k != l) {

          theta_block <- theta_nj_transpose[[l,k]]
          block_resid <- block_resid + (theta_block %*% Theta[-j,j][[l]] %*% SampleCov[[j,j]])

        }
      }

      block_resid <- block_resid + SampleCov[-j,j][[k]]
      fNorm <- norm(block_resid, type = "F")

      det(block_resid)

      if (fNorm <= gamma) {

        w_j <- matrix(0, nrow = M, ncol = M)

      } else {

        w_j <- matrix(fsolve(f = get_w_jk, block_resid = block_resid,  maxiter = 1000, x = diag(M),
                             Theta_inv = Theta_inv, SampleCov = SampleCov, gamma = gamma,
                             M = M, j = j, k = k)$x, nrow = M, ncol = M)
      }
    }


    for (block in 1:length(Theta[-j,j])) {

      Theta[-j,j][[block]] <- w_j[[block]]

    }

    error[iteration] <- norm(do.call(rbind, last_w) - do.call(rbind, Theta[-j,j]), type = "F")
    
    print(paste("Iteration ---", iteration, "--- Error ---", error[iteration]))
    # print(paste("iteration ---", iteration, "time:", Sys.time() - iter_time))
    iteration = iteration + 1

  }
  
  # print("W converged --- Leaving Algorithm 3")
  return(Theta)
  
}

##################################################################################
###    get wjk by solving system of equations                                  ###
###    (called by Algorithm_3)                                                 ###
##################################################################################
get_w_jk <- function(x, block_resid, Theta_inv, SampleCov, gamma, M = M, j = j, k = k) {
  
  # ans <- rep(0, M^2)
  w <- matrix(x, nrow = M, ncol = M)
  
  kron_prod <- kronecker(Theta_inv[-j,-j][[k,k]], SampleCov[[j,j]])

  return( c( (kron_prod %*% as.vector(w)) +
              as.vector(block_resid) + (gamma * (as.vector(w) / norm(w, type = "F"))) ) )
  
  # return(c(Theta_inv[-j,-j][[k,k]] %*% w %*% SampleCov[[j,j]] + block_resid + gamma * (w / norm(w, type = "F"))))
}


##################################################################################
###    get Uj for updating Algorithm 1 part (c)                                ###
###    (called by Algorithm_1)                                                 ###
##################################################################################
get_U <- function(Theta_inv, Theta, p, M, j) {
  
  inner_block <- list(matrix(0, nrow = M, ncol = M))
  
  
  Uj <- foreach(rowblock = 1:(p-1)) %dopar% {

    U_block <- matrix(0, nrow = M, ncol = M)

    foreach(colblock = 1:(p-1)) %do% {
      
      U_block <- U_block + (Theta_inv[-j,-j][[rowblock,colblock]] %*% Theta[-j,j][[colblock]])
    }
    
    Uj <- U_block
    }

  return(Uj)
}


##################################################################################
###    Update Sigma --- algorithm 1 part (c)                                   ###
###    (called by Algorithm_1)                                                 ###
##################################################################################
updateSigma <- function(Sigma, Theta_inv, sampleCov, Uj, p, j) {

  Uj_transpose <- lapply(Uj, t)

  # (1) update Sigma_jj
  Sigma[[j,j]] <- sampleCov[[j,j]]
  
  # (2) update small_sigma_j (multiply by -1 first or last????????????????????)
  for (block in 1:length(Uj)) {

    Sigma[-j,j][[block]] <- (-1*(Uj[[block]])) %*% sampleCov[[j,j]]
    
  }
  # (2.1) transpose of small_sigma_j (double check to make sure this calculation is correct)
  Sigma[j,-j] <- t(lapply(Sigma[-j,j], t))
  
 # (3) update Sigma 
  temp_sigma <- foreach(row_block = 1:(p-1)) %dopar% {

    foreach(col_block = 1:(p-1)) %do% {

      temp_sigma  <- (Theta_inv[-j,-j][[row_block, col_block]]
                      + Uj[[row_block]] %*% sampleCov[[j,j]] %*% Uj_transpose[[col_block]])

    }
  }
  
  # (3.1) update Sigma_j from list output of above parallel loops
  for (row_block in 1:length(temp_sigma)) {
    
    for (col_block in 1:length(temp_sigma[[row_block]])) {
      
      Sigma[-j,-j][[row_block, col_block]] <- temp_sigma[[row_block]][[col_block]]
      
    }
  }
  
  return(Sigma)
}


##################################################################################
###    Get Frobenius norm of current and previous precision matrix Theta       ###
###    (called by Algorithm_1)                                                 ###
##################################################################################
get_Norm <- function(Sigma, last_sigma, p) {
  
  unblock_sigma <- vector("list", p)
  unblock_last_sigma <- vector("list", p)
  
  for (row in 1:p) {
    
    unblock_sigma[[row]] <- do.call(cbind, Sigma[row,])
    unblock_last_sigma[[row]] <- do.call(cbind, last_sigma[row,])
    
  }
  
  unblock_sigma <- do.call(rbind, unblock_sigma)
  unblock_last_sigma <- do.call(rbind, unblock_last_sigma)
  
  fNorm <- norm((unblock_last_sigma - unblock_sigma), type = "F")
  return(fNorm)
}

##################################################################################
###    Get Edges from Precision Matrix                                         ###
###    (called by Algorithm_1)                                                 ###
###    Returns the final network structure                                     ###
##################################################################################
get_Precision <- function(Theta, SampleCov, Theta_inv, p, M) {
  start_time <- Sys.time()
  
  for (j in 1:p) {
    
    w_j_transpose <- t(matrix(lapply(Theta[-j,j], t), ncol = (p-1)))
    unblock_theta_nj <- vector("list", p-1)
    
    for (row in 1:(p-1)) {
      
      unblock_theta_nj[[row]] <- do.call(cbind, Theta_inv[-j,-j][row,])
      
    }
    
    unblock_theta_nj <- do.call(rbind, unblock_theta_nj)
    Theta[[j,j]] <- solve(SampleCov[[j,j]]) + (do.call(cbind, w_j_transpose) %*% unblock_theta_nj %*% do.call(rbind, Theta[-j,j]))
    print(Theta[[j,j]])
  }
  
  print(Sys.time() - start_time) 
  return(Theta)
  
}
  
  










  
  