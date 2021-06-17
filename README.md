# Introduction to BOSO


We present *BOSO*, an R package to perform feature selection in a linea regression problem. It implements a Bilevel Optimization Selector of Operators.

## Installation
BOSO can be installed from CRAN repository:

`install.packages("BOSO")`


## Introduction

The package package has been prepared to work like glmnet and lasso, presented 
in the *BestSubset* package.

``` r
library(BOSO)

## Load the data prepared for this test
load(system.file("data/high-5.rda", package = "BOSO"))

Xtr <- sim.xy$x
Ytr <- sim.xy$y
Xval <- sim.xy$xval
Yval <- sim.xy$yval


## Perform BOSO
time <- Sys.time()
obj <- BOSO(x = Xtr, y = Ytr,
            xval = Xval, yval = Yval,
            IC = 'eBIC',
            nlambda=100,
            intercept= 0,
            standardize = 0,
            Threads=4, timeLimit = 60, verbose = 3, 
            seed = 2021)
time <- as.numeric(Sys.time() - time)

```

`obj` is a BOSO object, which have the following associated functions: 

  - `coef(obj)` returns the coefficients (betas) of the linnear regression.  
  - `predict(obj, xnew)` returns the predicted outcome with a new X matrix.


``` r
betas <- coef(obj)
print(betas[betas!=0])

Ytr_predicted <- predict(obj, Xtr)
print(paste0("MSE for training set is ",  round(mean((Ytr_predicted-Ytr)^2),5)))

Yval_predicted <- predict(obj, Xval)
print(paste0("MSE for validation set is ", round(mean((Yval_predicted-Yval)^2),5)))
```








