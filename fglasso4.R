# fglasso4.R
# Perform fglasso simulations

###################################################################################
###         fglasso = ___Main___ Algorithm                                      ###
###                                                                             ###
### (1) 'iterbreak', 'epibreak', 'itertail' are essentially a tuple that        ###
### allows the algorithm to return a "second best" approxmination of the        ###
### precision matrix. When the 'epsilon' convergence criterion is not met       ###
### after at least 'iterbreak' iterations, but the previous 'itertail'          ###
### iteration errors are less than 'epibreak', the algorithm will return        ###
### the current precision matrix ,                                              ###
###                                                                             ###
### (2) 'iterlimit' is a stop early failure error. The algorithm will stop      ###
### running and return NaN and error information after 'iterlimit' iterations   ###
### with no 'epsilon' convergence criterion met.                                ###
###################################################################################
fglasso <- function(bmatrix, gamma, epsilon = 0.001, cores = 1, doprint = FALSE, 
                    iterbreak = 50, epibreak = (epsilon)*10, itertail = 5,
                    iterlimit = 100) {
  
  start_time <- as.numeric(Sys.time())
  
  require(nleqslv)
  
  if (cores == 0) {
    warning("Argument 'cores' cannot be zero")
  } else if (cores != 1) {
    require(doParallel)
    if (cores > detectCores()) {
      warning("Argument 'cores' must be no larger than the number of cores on the computers CPU")
    } else if (cores < 0){
      registerDoParallel(cores = detectCores() - abs(cores))
      print(paste("Number Cores in Use ---", getDoParWorkers()))
    } else {
      registerDoParallel(cores = cores)
      print(paste("Number Cores in Use ---", getDoParWorkers()))
    }
  }
  
  conv_fail <- FALSE
  
  p <- dim(bmatrix)[1]
  M <- unique(sapply(bmatrix, nrow))
  
  inner_block <- list(matrix(0, nrow = M, ncol = M))
  Sigma <- matrix(inner_block, nrow = p, ncol = p)
  Theta <- matrix(inner_block, nrow = p, ncol = p)
  diag(Sigma) <- list(diag(M))
  diag(Theta) <- list(diag(M))
  Theta_inv <- Theta
  
  norm_error <- c(epsilon+10^3)
  iteration <- 1
  
  while (tail(norm_error, n = 1) > epsilon) {
    
    iter_time <- as.numeric(Sys.time())
    last_sigma <- Sigma
    broken <- FALSE
    
    for (j in 1:p) {
      if (doprint == TRUE) {
        j_time <- as.numeric(Sys.time())
        print(paste("Function: ", j, " / ", p))
        print(paste("gamma:", gamma))
      }
      Theta_inv <- updateTheta(Theta = Theta, Sigma = Sigma, Theta_inv = Theta_inv, p = p, j = j, cores = cores)
      Theta <- Algorithm_3(Theta = Theta, Theta_inv = Theta_inv, SampleCov = bmatrix, 
                           gamma = gamma, epsilon = epsilon, p = p, M = M, j = j, 
                           cores = cores, doprint = doprint)
      if (length(Theta) == 0) {
        broken <- TRUE
        break
      }
      Uj <- get_U(Theta_inv = Theta_inv, Theta = Theta, p = p, M = M, j = j, cores = cores)
      Sigma <- updateSigma(Sigma = Sigma, Theta_inv = Theta_inv, sampleCov = bmatrix, 
                           Uj = Uj, p = p, j = j, cores = cores)
      if (doprint == TRUE) {
        print(paste("Function:", j, " time:", as.numeric(Sys.time()) - j_time))
      }  
    }
    if (broken == TRUE) {
      conv_fail <- TRUE
      break
    } else {
      fNorm <- get_Norm(Sigma = Sigma, last_sigma = last_sigma, p = p)
      norm_error[iteration] <- fNorm
      plot(seq(1:iteration), norm_error, xlab = "Iteration", ylab = "Error")
      if (doprint == TRUE) {
        print("--------------------------------------------------")
        print("--------------------------------------------------")
        print(paste("ITERATION:", iteration, " ERROR ---", norm_error[iteration]))
        print(paste("ITERATION:", iteration, " TIME ---", as.numeric(Sys.time()) - iter_time))
        print("--------------------------------------------------")
        print("--------------------------------------------------")
      }
      
      iteration <- iteration + 1
      
      if (iteration >= iterbreak & all(tail(norm_error, itertail) <= epibreak)) {
        PrecisionMatrix <- get_Precision(Theta = Theta, SampleCov = bmatrix, Theta_inv = Theta_inv, p = p, M = M)
        edges <- getEdges(PrecisionMatrix)
        return(list("edges" = edges, "PrecisionMatrix" = PrecisionMatrix, "Time" = as.numeric(Sys.time()) - start_time, 
                    "numIterations" = iteration, "normErrors" = norm_error, "finalNormError" = fNorm))
        
      } else if (iteration >= iterlimit) {
        print("Algorithm failed to converge: error 1")
        return(list("PrecisionMatrix" = NaN, "Time" = as.numeric(Sys.time()) - start_time, 
                    "numIterations" = iteration, "normErrors" = norm_error, "errorNum" = 1))
      }
    }
  }
  if (conv_fail == TRUE) {
    if (iteration == 1) {
      norm_error <- NaN
    }
    print("Algorithm failed to converge: error 2")
    return(list("PrecisionMatrix" = NaN, "Time" = as.numeric(Sys.time()) - start_time, 
                "numIterations" = iteration, "normErrors" = norm_error, "errorNum" = 2))
  } else {
    print(paste("Time to Convergence ---", as.numeric(Sys.time()) - start_time))
    PrecisionMatrix <- get_Precision(Theta = Theta, SampleCov = bmatrix, Theta_inv = Theta_inv, p = p, M = M)
    edges <- getEdges(PrecisionMatrix)
    return(list("edges" = edges, "PrecisionMatrix" = PrecisionMatrix, "Time" = as.numeric(Sys.time()) - start_time, 
                "numIterations" = iteration, "normErrors" = norm_error, "finalNormError" = fNorm))
  }
}

############################################################################
###           Step 2(a) of algorithm 1                                   ###
############################################################################
updateTheta <- function(Theta, Sigma, Theta_inv, p, j, cores) {
  
  SigSmallTranspose <- t(matrix(lapply(Sigma[-j,j], t), ncol = length(Sigma[-j,j])))
  
  if (cores != 1) {
    temp_theta <- foreach(row_block = 1:(p-1)) %dopar% {
      foreach(col_block = 1:(p-1)) %do% {
        temp_theta  <- (Sigma[-j,-j][[row_block, col_block]]
                        - Sigma[-j,j][[row_block]] %*% solve(Sigma[[j,j]]) %*% SigSmallTranspose[[col_block]])
      }
    }
    for (rowblock in 1:length(temp_theta)) {
      Theta_inv[-j,-j][rowblock,] <- temp_theta[[rowblock]]
    }
  } else {
    for (row_block in 1:(p-1)) {
      for (col_block in 1:(p-1)) {
        Theta_inv[-j,-j][[row_block, col_block]] <- (Sigma[-j,-j][[row_block, col_block]]
                                                     - Sigma[-j,j][[row_block]] %*% solve(Sigma[[j,j]]) %*% SigSmallTranspose[[col_block]])
      }
    }
  }
  return(Theta_inv)
}

##################################################################################
###    Algorithm 3: solve for w_j                                              ###
##################################################################################
Algorithm_3 <- function(Theta, Theta_inv, SampleCov, gamma, epsilon, p, M, j, cores, doprint) {
  
  iteration = 1
  error <- c(epsilon+10^3)
  broken <- FALSE
  
  while (tail(error, n = 1) > epsilon) {
    if (doprint == TRUE) {
      print(paste("w iteration ---", iteration))
    }
    last_w <- Theta[-j,j]
    theta_nj_transpose <- t(matrix(lapply(Theta_inv[-j,-j], t), nrow = nrow(Theta_inv[-j,-j]), ncol = ncol(Theta_inv[-j,-j])))
    
    if (cores != 1) {  ### PARALLEL ###
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
        if (fNorm <= gamma) {
          w_j <- matrix(0, nrow = M, ncol = M)
        } else {
          w_j <- matrix(nleqslv(x = diag(M), fn = get_w_jk,
                                jac = NULL, control=list(btol=.01, allowSingular = TRUE, maxit = 100),
                                Theta_inv = Theta_inv, SampleCov = SampleCov,
                                block_resid = block_resid, gamma = gamma,
                                M = M, j = j, k = k)$x, nrow = M, ncol = M)
        }
      }
      Theta[-j,j] <- w_j
      Theta[j,-j] <- lapply(w_j, t)
      
    } else {  ### NON PARALLEL ###
      
      w_j <- vector("list", (p-1))
      for(k in 1:(p-1)) {
        block_resid <- matrix(0, nrow = M, ncol = M)
        for(l in 1:(p-1)){
          if (k != l) {
            theta_block <- theta_nj_transpose[[l,k]]
            block_resid <- block_resid + (theta_block %*% Theta[-j,j][[l]] %*% SampleCov[[j,j]])
          }
        }
        block_resid <- block_resid + SampleCov[-j,j][[k]]
        fNorm <- norm(block_resid, type = "F")
        if (fNorm <= gamma) {
          w_j[[k]] <- matrix(0, nrow = M, ncol = M)
        } else {
          # w_j[[k]] <- matrix(nleqslv(x = diag(M), fn = get_w_jk,
          #                            jac = NULL, control=list(btol=.01, allowSingular = TRUE, maxit = 100),
          #                            Theta_inv = Theta_inv, SampleCov = SampleCov,
          #                            block_resid = block_resid, gamma = gamma,
          #                            M = M, j = j, k = k)$x, nrow = M, ncol = M)
          initial_guess <- matrix(solve(kronecker(Theta_inv[-j,-j][[k,k]], SampleCov[[j,j]]))%*%(-c(block_resid)), M, M)
          w_j[[k]] <- matrix(nleqslv(x = initial_guess,fn = get_w_jk, jac = NULL, 
                                     control=list(btol=.01, allowSingular = TRUE, maxit = 100),
                                     Theta_inv = Theta_inv, SampleCov = SampleCov,
                                     block_resid = block_resid, gamma = gamma,
                                     M = M, j = j, k = k)$x, nrow = M, ncol = M)
        }
      }
      Theta[-j,j] <- w_j
      Theta[j,-j] <- lapply(w_j, t)
    }
    error[iteration] <- norm(do.call(rbind, last_w) - do.call(rbind, Theta[-j,j]), type = "F")
    iteration = iteration + 1
    if (iteration > 100) {
      broken <- TRUE
      break
    }
  }
  if (broken == TRUE) {
    return(list())
  } else {
    return(Theta)
  }
}

##################################################################################
###    get wjk by solving system of equations                                  ###
###    (called by Algorithm_3)                                                 ###
##################################################################################
get_w_jk <- function(x, Theta_inv, SampleCov,
                     block_resid, gamma, M, j, k) {
  
  w <- matrix(x, nrow = M, ncol = M)
  kron_prod <- kronecker(Theta_inv[-j,-j][[k,k]], SampleCov[[j,j]])
  wj <- ((kron_prod %*% as.vector(w)) +
           as.vector(block_resid) + (gamma * (as.vector(w) / norm(w, type = "F"))))
  
  return(wj)
}

##################################################################################
###    get Uj for updating Algorithm 1 part (c)                                ###
###    (called by fglasso)                                                 ###
##################################################################################
get_U <- function(Theta_inv, Theta, p, M, j, cores) {
  
  if (cores != 1) {
    Uj <- foreach(rowblock = 1:(p-1)) %dopar% {
      U_block <- matrix(0, nrow = M, ncol = M)
      foreach(colblock = 1:(p-1)) %do% {
        U_block <- U_block + (Theta_inv[-j,-j][[rowblock,colblock]] %*% Theta[-j,j][[colblock]])
      }
      Uj <- U_block
    }
    
  } else {
    inner_block <- list(matrix(0, nrow = M, ncol = M))
    Uj <- rep(inner_block, length(Theta[-j,j]))
    for (rowblock in 1:length(Uj)) {
      for(colblock in 1:length(Uj)) {
        Uj[[colblock]] <- Uj[[colblock]] + (Theta_inv[-j,-j][[rowblock,colblock]] %*% Theta[-j,j][[colblock]])
      }
    }
  }
  return(Uj)
}

##################################################################################
###    Update Sigma --- algorithm 1 part (c)                                   ###
###    (called by fglasso)                                                 ###
##################################################################################
updateSigma <- function(Sigma, Theta_inv, sampleCov, Uj, p, j, cores) {
  
  Uj_transpose <- lapply(Uj, t)
  # (1) update Sigma_jj
  Sigma[[j,j]] <- sampleCov[[j,j]]
  # (2) update small_sigma_j 
  for (block in 1:length(Uj)) {
    Sigma[-j,j][[block]] <- -1*(Uj[[block]] %*% sampleCov[[j,j]])
  }
  # (2.1) transpose of small_sigma_j 
  Sigma[j,-j] <- t(lapply(Sigma[-j,j], t))
  # (3) update Sigma 
  if (cores != 1) {
    temp_sigma <- foreach(row_block = 1:(p-1)) %dopar% {
      foreach(col_block = 1:(p-1)) %do% {
        temp_sigma  <- (Theta_inv[-j,-j][[row_block, col_block]]
                        + Uj[[row_block]] %*% sampleCov[[j,j]] %*% Uj_transpose[[col_block]])
      }
    }
    # (3.1) update Sigma_j from list output of above parallel loops
    for (rowblock in 1:length(temp_sigma)) {
      Sigma[-j,-j][rowblock,] <- temp_sigma[[rowblock]]
    }
  } else {
    for (row_block in 1:(p-1)) {
      for (col_block in 1:(p-1)) {
        Sigma[-j,-j][[row_block, col_block]] <- (Theta_inv[-j,-j][[row_block, col_block]]
                                                 + Uj[[row_block]] %*% sampleCov[[j,j]] %*% Uj_transpose[[col_block]])
      }
    }
  }
  return(Sigma)
}

##################################################################################
###    Get Frobenius norm of current and previous precision matrix Theta       ###
###    (called by fglasso)                                                 ###
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
###             Get Final Precision Matrix                                     ###
##################################################################################
get_Precision <- function(Theta, SampleCov, Theta_inv, p, M) {
  
  for (j in 1:p) {
    w_j_transpose <- t(matrix(lapply(Theta[-j,j], t), ncol = (p-1)))
    unblock_theta_nj <- vector("list", p-1)
    for (row in 1:(p-1)) {
      unblock_theta_nj[[row]] <- do.call(cbind, Theta_inv[-j,-j][row,])
    }
    unblock_theta_nj <- do.call(rbind, unblock_theta_nj)
    Theta[[j,j]] <- solve(SampleCov[[j,j]]) + (do.call(cbind, w_j_transpose) %*% unblock_theta_nj %*% do.call(rbind, Theta[-j,j]))
  }
  return(Theta)
}

##################################################################################
###    Get estimated edges from Precision Matrix                               ###
##################################################################################
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