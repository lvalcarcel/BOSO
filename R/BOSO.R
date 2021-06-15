#' BOSO and associates functions
#'
#' Compute the BOSO for ust one block. This function calls \code{\link{cplexAPI}}
#' to solve the optimization problem
#' 
#'
#' @param x Input matrix, of dimension 'n' x 'p'. This is the data from the 
#' training partition. Its recommended to be class "matrix".
#' 
#' @param y Response variable for the training dataset. A matrix of one column 
#' or a vector, with 'n' elements.
#' 
#' @param xval Input matrix, of dimension 'n' x 'p'. This is the data from the 
#' validation partition. Its recommended to be class "matrix".
#' 
#' @param yval Response variable for the validation dataset. A matrix of one 
#' column or a vector, with 'n' elements.
#' 
#' @param intercept Boolean variable to indicate if intercept should be added 
#' or not. Default is false.
#' 
#' @param standardize Boolean variable to indicate if data should be scaled 
#' according to mean(x) mean(y) and sd(x) or not. Default is false.
#' 
#' @param dfmax Maximum number of variables to be included in the problem. The 
#' intercept is not included in this number. Default is min(p,n).
#' 
#' @param maxVarsBlock maximum number of variables in the block strategy.
#' 
#' @param metric information criteria to be used. Default is 'eBIC'.
#' 
#' @param nlambda The number of lambda values. Default is 100.
#' 
#' @param nlambda.blocks The number of lambda values in the block strategy part. 
#' Default is 10.
#' 
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of 
#' lambda.max, the (data derived) entry value.
#' 
#' @param lambda A user supplied lambda sequence. Typical usage is to have the 
#' program compute its own lambda sequence based on nlambda and 
#' lambda.min.ratio. Supplying a value of lambda overrides this. 
#' WARNING: use with care.
#'  
#' @param costErrorVal Cost of error of the validation set in the objective 
#' function. Default is 1. WARNING: use with care, changing this value changes 
#' the formulation presented in the main article.
#' 
#' @param costErrorTrain Cost of error of the training set in the objective 
#' function. Default is 0. WARNING: use with care, changing this value changes 
#' the formulation presented in the main article.
#' 
#' @param costVars Cost of new variables in the objective function. Default is 0. 
#' WARNING: use with care, changing this value changes the formulation 
#' presented in the main article.
#' 
#' @param Threads CPLEX parameter, number of cores that cplex is allowed to use.
#' Default is 0 (automatic).
#' 
#' @param timeLimit CPLEX parameter, time limit per problem provided to CPLEX.
#' Default is 1e75 (infinite time).
#' 
#' @param verbose print progress. Default is FALSE.
#' 
#' @param seed set seed for random number generator for the block strategy. 
#' Default is system default.
#' 
#' @param warmstart warmstart for CPLEX or use a different problem for each k. 
#' Default is False.
#' 
#' @param TH_metric is the ratio over one that the information criteria must 
#' increase to be STOP. Default is 1e-2.
#' 
#' @param indexSelected array of pre-selected variables. WARNING: debug feature.
#' 
#' @description Bonjour
#'
#' @author Luis V. Valcarcel
#' @export BOSO

BOSO = function(x, y, xval, yval, metric = 'eBIC',
                nlambda=100, nlambda.blocks = 10,
                lambda.min.ratio=ifelse(nrow(x)<ncol(x),0.01,0.0001),
                lambda=NULL, intercept=TRUE, standardize=TRUE,
                dfmax = NULL,
                maxVarsBlock = 10,
                costErrorVal = 1, costErrorTrain = 0, costVars = 0,
                Threads=0, timeLimit = 1e75, verbose = F,
                seed = NULL,
                warmstart = F,
                TH_metric = 1e-2,
                indexSelected = NULL)  {
  
  # Check for cplexAPI package
  if (!requireNamespace('cplexAPI', quietly = T)) {
    stop("Package cplexAPI not installed (required here)!")
  }
  
  # rm (list = ls())
  # setwd("D:/PhD/3 - Machine Learning MILP/LinearRegression-2021-02/R_package/BOSO R/BOSO")
  # source("R/BOSO_multiple_ColdStart.R")
  # source("R/BOSO_multiple_WarmStart.R")
  # source("R/utils.R")
  # sim.xy <- readRDS("D:/PhD/3 - Machine Learning MILP/Datasets/Hastie/RDSfiles/sim.xy.n500.p100.beta2.rho0.00.snr0.05.rds")
  # sim.xy <- sim.xy[[1]]
  # 
  # intercept <- F
  # standardize <- F
  # warmstart <- T
  # 
  # if (intercept|standardize) {
  #   sim.xy$betas <- c(1, sim.xy$beta)
  # } else
  #   sim.xy$betas <- sim.xy$beta
  # 
  # 
  # x <- sim.xy$x
  # y <- sim.xy$y
  # xval <- sim.xy$xval
  # yval <- sim.xy$yval
  # 
  # metric = 'eBIC'
  # nlambda=100
  # nlambda.blocks = 10
  # lambda.min.ratio=ifelse(nrow(x)<ncol(x),0.01,0.0001)
  # lambda=NULL
  # # intercept=T
  # # standardize=F
  # dfmax = NULL
  # maxVarsBlock = 5
  # costErrorVal = 1
  # costErrorTrain = 0
  # costVars = 0
  # Threads=0
  # timeLimit = 1e75
  # verbose = 5
  # seed = NULL
  # TH_metric <- 1e-3
  
  # Set up data ####
  x = as.matrix(x)
  y = as.numeric(y)
  xval = as.matrix(xval)
  yval = as.numeric(yval)
  n = nrow(x)
  nval = nrow(xval)
  p = ncol(x)
  
  # standarze?
  if (standardize) {
    # standardize using xtrain and xval
    # obj = standardize(rbind(x, xval), c(y, yval), intercept=T, normalize = T)
    obj = standardize(x, y, intercept=T, normalize = T)
    x = obj$x
    y = obj$y
    mx = obj$mx
    my = obj$my
    sx = obj$sx
    # obj = standardize(x, y, mx=mx, my=my, sx=sx)
    # x = obj$x
    # y = obj$y
    obj = standardize(xval, yval, mx=mx, my=my, sx=sx)
    xval = obj$x
    yval = obj$y
    # add for later on
    intercept = F  # once data is scaled, no need for beta0
    mx = c(0,mx) 
    sx = c(1,sx)
  } else {
    mx = rep(0,p + 1)
    my = 0
    sx = rep(1,p + 1)
  }
  
  # Set dfmax manually if NULL
  if (is.null(dfmax)){dfmax = p} 
  
  # Generate the lambda array if necessary
  if (is.null(lambda)){
    # lambda_max <- norm(t(x)%*%y, "I")/nrow(x)
    lambda_max <- max( abs(t(y - mean(y)*(1-mean(y))) %*% x ) )/ n #lasso
    lambda_min <- lambda_max * lambda.min.ratio
    lambda <- exp(seq(log(lambda_max*1e3), log(lambda_min), length.out = nlambda)) # lambda_max from ridge, lambda_min from lasso
    lambda.blocks <- exp(seq(log(lambda_max*1e3), log(lambda_min), length.out = nlambda.blocks))
  } else {
    # Reset nlambda if a specific lambda sequence is passed
    nlambda <- length(lambda)
    nlambda.blocks <- length(lambda.blocks)
    lambda.blocks <- lambda.blocks
  }
  
  # if eBIC, set BIC for blocks
  metric.block <- ifelse(metric=="eBIC", "BIC", metric)
  
  
  ## Separate data in packages of with a maximum size and perform block strategy ####
  data.raw <- list(x = x, y = y, xval = xval, yval = yval)
  
  if (is.null(indexSelected)){
    indexSelected <- 1:p # index of variables
  }
  
  # set seed for random number generator
  if (!is.null(seed)){ set.seed(seed) }
  
  # auxiliary variables
  FinishProblem <- F
  ContinueBlocks <- T
  numIter <- 0
  indexSelectedIter <- list()
  
  
  # maxVarsBlock = 10
  # indexSelected <- 1:85
  
  while (ContinueBlocks & length(indexSelected)>maxVarsBlock*1.5){
    numIter <- numIter+1
    
    #separate into blocks
    idx <- sample(indexSelected)
    k <- ceiling(length(idx)/maxVarsBlock)*maxVarsBlock 
    if(k!=length(idx)) {idx[(length(idx)+1):k] <- NA}
    idx <- matrix(idx, nrow = maxVarsBlock, byrow = T)
    idx <- lapply(lapply(apply(idx,2,as.list),unlist),function(x){x[!is.na(x)]})  # transform each column in a list
    
    indexSelectedIter[[numIter]] <- list()
    time_block <- rep(NA, length(idx))
    
    if (verbose) { cat(paste0('Iteration ',numIter,'\t number of variables: ',length(indexSelected),' \t block size: ',maxVarsBlock,'\n')) }
    
    # solve each of the blocks
    for (i in 1:length(idx)){
      # i = 1
      p.block <- length(idx[[i]])
      numVarArray <- seq(0,p.block)

      # subset of variables from the block strategy
      x = data.raw$x[,idx[[i]]]
      xval = data.raw$xval[,idx[[i]]]
      
      time_block[i] <- Sys.time()
      
      if (verbose>1) {  cat(paste0('Iteration ',numIter,'\t number of variables: ', length(indexSelected), '\t subproblem ', round(i*100/length(idx)), '%\n' )) }
      
      if (warmstart) {
        # cat("\nWARMSTART\n")
        result.block <- BOSO.multiple.warmstart(x = x, y = y, xval = xval, yval = yval,
                                                lambda=lambda.blocks, 
                                                intercept=intercept, 
                                                standardize=F,
                                                dfmin = 0, dfmax = p.block,
                                                costErrorVal = costErrorVal, costErrorTrain = costErrorVal, 
                                                costVars = costVars,
                                                Threads=Threads, timeLimit = timeLimit, verbose = max(verbose-2,0),
                                                metric = metric.block, n.metric = n+nval, p.metric = p,
                                                TH_metric = TH_metric)
      } else {
        # cat("\ncold start\n")
        result.block <- BOSO.multiple.coldstart(x = x, y = y, xval = xval, yval = yval,
                                                lambda=lambda.blocks, 
                                                intercept=intercept, 
                                                standardize=F,
                                                dfmin = 0, dfmax = p.block,
                                                costErrorVal = costErrorVal, costErrorTrain = costErrorVal, 
                                                costVars = costVars,
                                                Threads=Threads, timeLimit = timeLimit, verbose = max(verbose-2,0),
                                                metric = metric.block, n.metric =  n+nval, p.metric = p,
                                                TH_metric = TH_metric)
      }
      
      if (nrow(result.block$betas)>length(idx[[i]])) {result.block$betas <- result.block$betas[-1,]}
      indexSelectedIter[[numIter]][[i]] <- idx[[i]][result.block$betas[,which.min(result.block$score)]!=0]
      
      time_block[i] <- as.numeric(Sys.time() - time_block[i])
      
      if (verbose>1) { cat(paste0('Iteration ',numIter,'\t number of variables: ',length(idx[i]),
                                  '\t block size: ',maxVarsBlock,'\tSelected = ',length(indexSelectedIter[[numIter]][[i]]),
                                  '\tElapsed time= ',round(time_block[i],3),'\n'))}
      
    } 
    indexSelectedIter[[numIter]] <- sort(unlist(indexSelectedIter[[numIter]]))
    indexSelected <- indexSelectedIter[[numIter]]
    
    if (verbose>=3){ print(indexSelected) }
    
    # check if we need to stay in the loop
    if (length(indexSelected) <= maxVarsBlock*1.5){
      ContinueBlocks = F
    } else if (numIter>2) {
      if(length(indexSelectedIter[[numIter]])<length(indexSelectedIter[[numIter-1]])){
        ContinueBlocks = T
      } else if (length(intersect(indexSelected,indexSelectedIter[[numIter-1]]))==length(indexSelected) &
                 length(intersect(indexSelected,indexSelectedIter[[numIter-2]]))==length(indexSelected)){
        ContinueBlocks = F
      }
    }
  }
  
  
  ## if there is no variable after block strategy (eg: low SNR), exit the problem ####
  
  if (length(indexSelected)==0){
    result <- list(x = x,
                   y = y,
                   xval = xval,
                   yval = yval, 
                   metric = metric,
                   nlambda = nlambda,
                   lambda = lambda,
                   intercept = intercept,
                   standardize = standardize,
                   mx = mx, my = my, sx = sx,
                   dfmax = dfmax,
                   lambda.selected = 0,
                   p = p, n = n, nval = nval)
    
    result$betas <- c(mean(y), rep(0, p))
    if (!intercept) {result$betas <- result$betas[-1]}
    result$betas <- matrix(result$betas, ncol = 1)
    
    if (intercept){
      result$errorTrain = y-mean(y)
      result$errorVal = yval-mean(y)
    } else {
      result$errorTrain = y 
      result$errorVal = yval
    }
    
    class(result) = "BOSO"
    object <- result
    return(result)
  }
  
  
  
  ### Final problem ####
  
  
  dfmax.raw <- dfmax
  p.final <- length(indexSelected)
  dfmax <-  ifelse(p.final > dfmax , dfmax, p.final)
  numVarArray <- seq(0,dfmax)
  # if (!intercept) {numVarArray <- numVarArray[-1]}
  
  # subset of variables from the block strategy
  x = data.raw$x[,indexSelected]
  xval = data.raw$xval[,indexSelected]
  
  if (verbose) {  cat(paste0('Final problem:\t number of variables: ', dim(x)[2], ' \t \n' )) }
  
  if (warmstart) {
    # cat("\nWARMSTART\n")
    result.final <- BOSO.multiple.warmstart(x = x, y = y, xval = xval, yval = yval,
                                            lambda=lambda, 
                                            intercept=intercept, 
                                            standardize=F,
                                            dfmin = 0, dfmax = dfmax,
                                            costErrorVal = costErrorVal, costErrorTrain = costErrorVal, 
                                            costVars = costVars,
                                            Threads=Threads, timeLimit = timeLimit, verbose = max(verbose-1,0),
                                            metric = metric, n.metric = n+nval, p.metric = p,
                                            TH_metric = TH_metric)
  } else {
    # cat("\ncold start\n")
    result.final <- BOSO.multiple.coldstart(x = x, y = y, xval = xval, yval = yval,
                                            lambda=lambda, 
                                            intercept=intercept, 
                                            standardize=F,
                                            dfmin = 0, dfmax = dfmax,
                                            costErrorVal = costErrorVal, costErrorTrain = costErrorVal, 
                                            costVars = costVars,
                                            Threads=Threads, timeLimit = timeLimit, verbose = max(verbose-1,0),
                                            metric = metric, n.metric =  n+nval, p.metric = p,
                                            TH_metric = TH_metric)
  }
  
  
  # object <- obj
  
  idx = which.min(result.final$score)
  
  # Append a few things to the returned object
  result <- list(x = x,
                 y = y,
                 xval = xval,
                 yval = yval, 
                 metric = metric,
                 nlambda = nlambda,
                 lambda = lambda,
                 intercept = intercept,
                 standardize = standardize,
                 mx = mx, my = my, sx = sx,
                 dfmax = dfmax,
                 result.final = result.final, 
                 errorTrain = result.final$errorTrain[,idx], 
                 errorVal = result.final$errorVal[,idx],
                 lambda.selected = result.final$lambda.selected[idx],
                 p = p, n = n, nval = nval,
                 blockStrategy = indexSelectedIter)
  
  
  result$betas <- rep(0, p + 1)
  result$betas[c(1,indexSelected+1)] <- result.final$betas[,idx]
  
  result$betas <- matrix(result$betas, ncol = 1)
  
  class(result) = "BOSO"
  
  # object <- result
  # cbind(sim.xy$betas, coef.BOSO(result))
  # table(real = sim.xy$betas!=0, BOSO = coef.BOSO(result)!=0)
  # c(object$lambda.selected, which(lambda==object$lambda.selected))
  # 
  # c(result.final$time, sum(result.final$time[result.final$time<1e5]))
  
  return(result)
}





