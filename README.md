# fglasso
Gaussian graphical model via functional graphical lasso for network analysis.

See 'run_fglasso.R" for demonstration. 

"fglasso4.R" - returns conditional dependencies of nodes/variables in (high dimensional) network data.

  - fglasso: performs fglasso algorithm on a blocked covariance matrix. 
  
                ARGUMENTS
                bmatrix = blocked covariance matrix in which each row/column represents a covaraince component 
                          of two nodes/variables.
                gamma = penalizing parameter to reduce covariance estimates to zero; a numeric. 
                epsilon = convergence norm error criterion; frobenus norm of current - previous precision matrix must be <
                          epsilon to converge. 
                cores = number of cpu cores to allocate for computations. NOTE: multiprocessing may slow down computation time 
                          for small dimensional problems. High dimensional problems display huge increases in processing speed 
                          via multiprocessing. 
                doprint = (boolean: TRUE/FALSE); print algorithm updates during computations. 
                iterbreak = max number of iterations to compute before returning suboptimal precision matrix (tuple with  
                          epibreak, itertail). 
                epibreak = secondary error criterion (numeric larger than epsilon) allowing for suboptimal convergence 
                          following 'iterbreak' iterations with error <= 'epibreak'. 
                itertail = number of consecutive iterations (following 'iterbreak' iterations) reaching 'epibreak' error 
                          criterion necessary before returning suboptimal precision matrix. 
                iterlimit = max number of iterations to compute before breaking the program and returning NaN. 
                

"jFoldPCA.R" - creates functional data objects and corresponding covariance matrix of fdobjs. 

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
