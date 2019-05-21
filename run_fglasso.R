# Scope: apply fglasso
# load 'sim_data.R' to get data
# load 'jFoldPCA.R' and run specified functions fit fdobjs and obtain covariance matrix
# load 'fglasso.R' and run Algorithm_1 to obtain functional graphical model 

#-----------------------------------------------------------#
#     1. load/simulate data.                                #
#-----------------------------------------------------------#
source("sim_data.R")

# simulated data: p x k x n array of simulated data.
simdat <- sim_data(n = 100, p = 100, prec_model = 1)

# turn 3D simulated data into list of 2D
simList <- vector("list", dim(simdat)[1])
for (j in 1:length(simList)) {
  simList[[j]] <- simdat[j,,]
}

#-----------------------------------------------------------#
#     2. determine best number of basis functions L.        #
#-----------------------------------------------------------#
source("jFoldPCA.R")

simFit <- fit_Basis(dataList = simList, nbasistry = 5:20)
plot(simFit$gcvs)

best_basis = simFit$best.basis
fdSim <- fit_Basis(dataList = simList, nbasistry = best_basis)$fdobjList

#-----------------------------------------------------------#
#     3. determine best number of harmonics by              #
#         performing J-fold cross-validation.               #
#-----------------------------------------------------------#
source("jFoldPCA.R")

cvpca <- jFoldPCA(dataList = simList, nbasis = best_basis, 
                  nfolds = 4, nharmrange = 1:best_basis)
plot(cvpca$nharm_sse)
best_nharm <- cvpca$best_nharm

#-----------------------------------------------------------#
#     3. Get sample covariance matrix of eigenscores        #
#         (using best_nharm)                                #
#         block the covariance matrix                       #
#-----------------------------------------------------------#
source("jFoldPCA.R")

# centerfns = TRUE or FALSE?

sample_cov <- pc_cov(fdobjList = fdSim, nharm = best_nharm, block = TRUE)

#-----------------------------------------------------------------------------------------------#
#     5. Begin fglasso algorithm                                                                #
#-----------------------------------------------------------------------------------------------#
source("fglasso4.R")

# fgm <- fglasso(bmatrix = sample_cov, gamma = 90, cores = 3, doprint = TRUE, iterbreak = 20, iterlimit = 25)

fgm_core1 <- fglasso(bmatrix = sample_cov, gamma = 110, cores = 1, doprint = TRUE, iterbreak = 20, iterlimit = 25) # failed at j=8 (g=110,p=100,n=100)
fgm_core3 <- fglasso(bmatrix = sample_cov, gamma = 110, cores = 3, doprint = TRUE, iterbreak = 20, iterlimit = 25)


# unblock_1 <- vector("list", 50)
unblock_3 <- vector("list", 50)
for (row in 1:50) {
  
  # unblock_1[[row]] <- do.call(cbind, fgm_core1$PrecisionMatrix[row,])
  unblock_3[[row]] <- do.call(cbind, fgm_core3$PrecisionMatrix[row,])
  
}
# unblock_1 <- do.call(rbind, unblock_1)
unblock_3 <- do.call(rbind, unblock_3)

# norm((unblock_1 - unblock_3), type = "F")

simulated <- get_precision(modelnum = 1, p = 100)
block_simulated <- block_matrix(inputmatrix = simulated, outersize = 50, innersize = best_nharm)

# estimated_1 <- getEdges(blocked_matrix = fgm_core1$PrecisionMatrix)
estimated_3 <- getEdges(blocked_matrix = fgm_core3$PrecisionMatrix)
simulated_edges <- getEdges(blocked_matrix = block_simulated)

# errRates1 <- getErrors(estimated_1, simulated_edges)
errRates3 <- getErrors(estimated_1, simulated_edges)

"Function:  15  /  50"

 # nlesqlv & pracma both fail oat J = 39 for prec_mod = 1 on fglasso.3. (for parallel) (non-parallel fails at j = 4)

# With gamma = 1:
# Warning messages:
#   1: In inv(A0) : Matrix appears to be singular.






#-----------------------------------------------------------#
#     5. Run fglasso with different gammas                  #
#-----------------------------------------------------------#
source("fglasso2.R")

# 14: In gamma_fgm[gamma] <- Algorithm_1(bmatrix = sample_cov,  ... :
                                         # number of items to replace is not a multiple of replacement length
gammas <- 10^seq(10, -2, length = 100)
gamma_fgm <- vector("list", length(gammas))
for (g in 1:length(gammas)) {
  # gamma <- gammas[g]
  print(paste("gamma ---", g, "/", length(gammas)))
  gamma_fgm[[g]] <- Algorithm_1(bmatrix = sample_cov, gamma = gammas[g])
}






