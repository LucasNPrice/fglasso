# fglasso
functional gaussian graphical model via functional graphical lasso for network analysis 

Work still in progress. 

"fglasso2.R" - fglasso to return conditional dependencies of node/variables in data; a network analysis algorithm. 

  - Algorithm_1: essentially __main__ function to perform the entire fglasso algorithm. 
  
                ARGUMENTS
                bmatrix = blocked covariance matrix in which each row/column represents a covaraince component 
                          of two nodes/variables.
                gamma = penalizing parameter to reduce covariance estimates to zero; a numeric. 
                

"jFoldPCA.R" to create functional objects. 

  - fit_Basis: performs cross validation to determine the number of basis functions to best represent data. 
  
               ARGUMENTS
               dataList = list of data; each element in list is an observation. 
               nbasistry = number of basis splines to attempt to fit functional objects 
               
               RETURNS
               best.basis = number of basis functions found to result in best fit of the data 
               gcvs = generalized cross-validated error for each number of basis functions fitted 
                      in the cross-validation process
               fdobjList = list of functional data objects representing the input dataList ; only 
                            returned if nbasistry = 1. 
              
  - jFoldPCA: performs j-fold cross validation principle components analysis to select optimal number
              of principle components to represenet data; creates fdobjs and performs j-fold cv pca on the fdobjs. 
              
              ARGUMENTS
              dataList = original data list to assess error of fitted values 
              nbasis = number of basis splines to be used to fit functional objects 
              nfolds = number of folds/groups for data splitting in cross-validation
              nharmrange = range of principle compnents/harmonics to try to fit. 
  
              RETURNS
              nharm_sse = list of sum squared error for each number of harmonics fit. 
              best_nharm = number of harmonic functions found to result in lowest cv error
              
  - pc_cov: fit functional objects and returns the covaraince matrix of the fitted objects
  
              ARGUMENTS
              fdobjList = list of functional data objects 
              nharm = number of principle compnents to retain and fit the fdobjs 
              block = whether or not to return a blocked matrix; default is FALSE
              
  - block_matrix: creates a blocked matrix 
  
              ARGUMENTS
              inputmatrix = input matrix which will become a blocked matrix 
              outersize = col/row length of full blocked matrix
              innersize = col/row length of each matrix within full blocked matrix 
              
              REUTNS
              blocked_matrix = original input matrix now in a blocked structure 