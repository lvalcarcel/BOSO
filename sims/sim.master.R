#' Master function for running simulations.
#'
#' Run a set of simulations with the specified configuration.
#'
#' @param n,p The number of training observations, and the number of predictors.
#' @param nval The number of validation observations.
#' @param reg.funs This is a list of functions, representing the regression
#'   procedures to be used (evaluated) in the simulation. Each element of the
#'   list must be a function that takes x, y (the training predictor matrix and
#'   response vector) as its only two (mandatory) arguments, and must return an
#'   object with associated coef and predict methods. The coef method must take
#'   obj (the returned object) and return a matrix of coefficients, with one
#'   column per tuning parameter value inherent to the regression method. The
#'   predict method must take obj, newx (the returned object and a new predictor
#'   matrix) and return a matrix of predictions, again with one column per
#'   tuning parameter value inherent to the regression method.
#' @param Number of repetitions of which to average the results. Default is 50.
#' @param seed Seed to be set for the overall random number generation, i.e.,
#'   set before repetitions are begun (for reproducibility of the simulation
#'   results). Default is NULL, which effectively sets no seed.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param file,file.rep Name of a file to which simulation results are saved
#'   (using saveRDS), and a number of repetitions after which intermediate
#'   results are saved. Setting file to NULL is interpreted to mean that no
#'   simulations results should be saved; setting file.rep to 0 is interpreted
#'   to mean that simulations results should be saved at the very end, i.e., no
#'   intermediate saving. Defaults are NULL and 5, respectively.
#' @param rho,s,beta.type,snr. Arguments to pass to \code{\link{sim.xy}}; see
#'   the latter's help file for details.
#'
#' @return A list with components err.train, err.val, err.test, prop, risk, nzs,
#'   fpos, fneg, F1, opt for the training error, validation error, test error,
#'   test proportion of variance explained, risk, number of selected nonzero
#'   coefficients, number of false positives, number of false negatives, F1
#'   measure, and relative optimism (difference in test error and training
#'   error, divided by training error), respectively.  These are each lists of
#'   length N, where N is the number of regression methods under consideration
#'   (the length of reg.funs). The ith element of each list is then a matrix of
#'   dimension nrep x m, where m is the number of tuning parameters inherent to
#'   the ith method.
#'
#' @seealso \code{\link{sim.xy}}
#' @author Trevor Hastie, Robert Tibshirani, Ryan Tibshirani
#' @references The structure of this simulation code based on that from the
#'   \code{conformalInference} package.
#' @example examples/ex.sim.master.R
#' @export sim.master

sim.master = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                      file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  
  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)
  
  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = 
    opt = runtime = vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) = names(F1) = 
    names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA,nrep)
  
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xval, yval
    xy.obj = sim.xy(n,p,nval,rho,s,beta.type,snr)
    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y,xy.obj$xval,xy.obj$yval)
        })[1]
        
        # Grab the estimated coefficients, and the predicted values on the
        # training and validation sets
        betahat = as.matrix(coef(reg.obj))
        m = ncol(betahat); nc = nrow(betahat)
        
        # Check for intercept
        if (nc == p+1) {
          intercept = TRUE
          betahat0 = betahat[1,]
          betahat = betahat[-1,]
        }
        else intercept = FALSE
        
        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =
            risk[[j]] = nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }
        
        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        delta = betahat - xy.obj$beta
        risk[[j]][i,] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i,] = risk[[j]][i,] + betahat0^2
        err.test[[j]][i,] = risk[[j]][i,] + xy.obj$sigma^2
        prop[[j]][i,] = 1 - err.test[[j]][i,] / err.null[i]
        nzs[[j]][i,] = colSums(betahat!=0)
        tpos = colSums((betahat!=0)*(xy.obj$beta!=0))
        fpos[[j]][i,] = nzs[[j]][i,]-tpos
        fneg[[j]][i,] = colSums((betahat==0)*(xy.obj$beta!=0))
        F1[[j]][i,] = 2*tpos/(2*tpos+fpos[[j]][i,]+fneg[[j]][i,])
        opt[[j]][i,] = (err.test[[j]][i,] - err.train[[j]][i,]) /
          err.train[[j]][i,]
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }
    
    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(bestsubset:::enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,fneg,F1,opt,sigma,runtime),file=file)
    }
  }
  
  # Save results now (in case of an error that might occur below)
  out = bestsubset:::enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,nzs,fpos,
               fneg,F1,opt,sigma,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  out = bestsubset:::choose.tuning.params(out)
  
  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}
