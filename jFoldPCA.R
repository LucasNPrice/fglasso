# jFoldPCA.R
# Create functional data objects and perform J-fold cross validation with PCA
# returns best number of basis functions best_basis
# return best number of pc functions/harmonics best_nharm

############################################################################
###                     Tuning for best nbasis                           ###
### nbasistry = sequence/vector of values to try for number basis        ###
### dataList = list of matrices; each matrix gets a separate fdobj       ###
############################################################################
fit_Basis <- function(dataList, nbasistry, doprint = FALSE) {
  
  require(fda)
  
  if ((length(unique(sapply(dataList, nrow))) > 1) == TRUE) {
    warning("At least one list is of different domain than others!")
  } else {
    domain <- c(1, unique(sapply(dataList, nrow)))
  }
  
  fdobjList <- vector("list", length(dataList))
  sum_gcv <- rep(0, length(nbasistry))

  for (l in 1:length(nbasistry)) {
    nb <- nbasistry[l]
    if (doprint == TRUE) {
      print(paste("Computing basis number: ", nb))
    }
    spline.basis = create.bspline.basis(rangeval = domain, nbasis = nb)
    for (n in 1:length(dataList)) {
      fdobjList[[n]] <- smooth.basis(argvals = min(domain):max(domain), y = dataList[[n]], fdParobj = spline.basis)
      sum_gcv[l] = sum_gcv[l] + sum(fdobjList[[n]]$gcv, na.rm = T)
    }
  }
  if (length(nbasistry) > 1) {
    nb <- nbasistry[which.min(sum_gcv)]
    return(list("best.basis" = nb, "gcvs" = sum_gcv))
  } else {
    return(list("fdobjList" = fdobjList, "nbasis" = nbasistry, "gcvs" = sum_gcv))
  }
}

#########################################################################################
###             J-Fold Cross Validation with PCA                                     ####
#########################################################################################

#---------------------------------------------------------------------------------------#
#           Step 1: cross validation of PCA error (called from jFoldPCA)                #
#           fdobjList = list of objects of class fdobj.                                 #
#           nharmrange = range of number of harmonics to be used.                       #
#           train_domain = training domain.                                             #
#           test_domain = testing domain (not to be confused with actual data values)   #
#---------------------------------------------------------------------------------------#
cv_fdPCA <- function(fdobjList, nharmrange, dataList, train_domain, test_domain, jfold = NULL, doprint = TRUE) {

  sse <- rep(0, length(nharmrange))
  for (k in 1:length(nharmrange)) {
    nharm = nharmrange[k]
    if (is.null(jfold)) {
      jfold = NULL
      print(paste("Computing harmonic number: ", nharm))
    } else {
      if (doprint == TRUE) {
        print(paste("Computing --- Fold: ", jfold,  " --- Harmonic: ", nharm))
      }
    }
    for (n in 1:length(fdobjList)) {
      train_pca <- pca.fd(fdobjList[[n]]$fd, nharm = nharm)
      g_hat <- eval.fd(test_domain, train_pca$harmonics) %*% t(train_pca$scores)
      error_diff <- g_hat - dataList[[n]][test_domain,]
      squared_error <- sum(sapply(error_diff, function(x) x^2))
      sse[k] <- sse[k] + squared_error
    }
  }
  return(sse)
}

#-----------------------------------------------------------------------------#
#         Step 2: J-Fold Cross Validation of PCA error                        #
#                 creates fdobjs and performs j-fold cv pca on them           #
#                 (calls cv_fdPCA() to perform J-Fold cv)                     #
#-----------------------------------------------------------------------------#
jFoldPCA <- function(dataList, nbasis, nfolds, nharmrange, doprint = TRUE) {
  
  if ((length(unique(sapply(dataList, nrow))) > 1) == TRUE) {
    warning("At least one list is of different domain than others!")
  } else {
    full_domain <- c(1, unique(sapply(dataList, nrow)))
  }
  spline.basis = create.bspline.basis(rangeval = full_domain, nbasis = nbasis)
  fdobjList <- vector("list", length(dataList))
  cvErrors <- rep(0, length(nharmrange))
  
  fold_size <- floor(max(full_domain) / nfolds)
  shuffle <- sample(min(full_domain):max(full_domain), max(full_domain))
  test_iter <- 1:fold_size
  
  for (fold in 1:nfolds) {
    vector("list", length(dataList))
    for (n in 1:length(dataList)) {
      train <- shuffle[-test_iter]
      test <- shuffle[test_iter]
      trainData <- dataList[[n]][train,]
      fdobjList[[n]] <- smooth.basis(argvals = train, y = trainData, fdParobj = spline.basis)
    }
    J_errors <- cv_fdPCA(fdobjList = fdobjList, nharmrange = nharmrange, 
                         dataList = dataList, train_domain = train, test_domain = test, jfold = fold, doprint = doprint)
    cvErrors <- cvErrors + J_errors
    test_iter <- test_iter + fold_size
  }
  return(list("nharm_sse" = cvErrors, "best_nharm" = nharmrange[which.min(cvErrors)]))
}

#------------------------------------------------------------------------------------#
#         Step 3: Re-run PCA on all data using best nharm from j-fold cv             #
#                 Get sample covariance/correlation matrix of the PC scores          #
#------------------------------------------------------------------------------------#

pc_cor <- function(fdobjList, nharm, method = "cor", block = FALSE) {
  
  fd_pc <- vector("list", length(fdobjList))
  
  for (i in 1:length(fd_pc)) {
    fd_pc[[i]] <- pca.fd(fdobjList[[i]]$fd, nharm = nharm, centerfns = TRUE)$scores
    dim(fd_pc[[i]]) # = eigenscores of jth function for all i 1:n samples; dim = n x nharm
  }
  eigenscores <- do.call(cbind, fd_pc) # = eigenscore matrix for all functions; dim = n x (p x nharm)
  
  if (method == "cor") {
    cor_matrix <- cor(eigenscores)
    if (block == FALSE) {
      return(cor_matrix)
    } else {
      source("block_matrix.R")
      cor_matrix <- block_matrix(inputmatrix = cor_matrix, outersize = length(fdobjList), innersize = nharm)
      return(cor_matrix)
    }
  } else if (method == "cov") {
    cov_matrix <- cov(eigenscores)
    if (block == FALSE) {
      return(cov_matrix)
    } else {
      source("block_matrix.R")
      cov_matrix <- block_matrix(inputmatrix = cov_matrix, outersize = length(fdobjList), innersize = nharm)
      return(cov_matrix)
    }
  } else {
    warning("method must be either 'cor' or 'cov'.")
  }
}

#------------------------------------------------------------------------------------#
#         pc_cov() is just a wrapper for pc_cor()                                    #
#------------------------------------------------------------------------------------#
pc_cov <- function(fdobjList, nharm, method = "cov", block = FALSE) {
  return(pc_cor(fdobjList = fdobjList, nharm = nharm, method = method, block = block))
}

#------------------------------------------------------------------------------------#
#         getEigenscores() simply calculates and returns eigenscores                 #
#------------------------------------------------------------------------------------#
getEigenscores <- function(fdobjList, nharm) {
  
  fd_pc <- vector("list", length(fdobjList))
  
  for (i in 1:length(fd_pc)) {
    fd_pc[[i]] <- pca.fd(fdobjList[[i]]$fd, nharm = nharm, centerfns = TRUE)$scores
  }
  eigenscores <- do.call(cbind, fd_pc)
  return(eigenscores)
}