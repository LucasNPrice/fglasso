# Scope: apply fglasso
# load 'sim_data.R' to get data
# load 'jFoldPCA.R' and run specified functions fit fdobjs and obtain covariance matrix
# load 'fglasso.R' and run Algorithm_1 to obtain functional graphical model 

#-----------------------------------------------------------#
#     1. load/simulate data.                                #
#-----------------------------------------------------------#
source("sim_data.R")

# simulated data: p x k x n array of simulated data.
simdat <- sim_data(n = 100, p = 50, prec_model = 1)

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

#-----------------------------------------------------------#
#     5. Begin fglasso algorithms                           #
#-----------------------------------------------------------#
source("fglasso2.R")

fgm <- Algorithm_1(bmatrix = sample_cov, gamma = 100)

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


# how to set error criterion ---- what is a "good" difference between matrices and what is a "good" number of 
# iterations meeting error criterion?

# parallelize algorithms
# bootstrap for optimal lambda 




