# Introduction to BOSO


We present *Rediscover*, an R package to identify mutually exclusive genomic events. It reimplements a privious R package (*Discover*) whose main contribution is a statistical analysis based on the Poisson-Binomial distribution that takes into account that some samples are more mutated than others. *Rediscover* is much faster than the <a href="https://github.com/NKI-CCB/DISCOVER" target="_blank"> discover </a> implementation.

## Installation
BOSO can be installed from CRAN repository:

`install.packages("BOSO")`


## Introduction

The package library has two main parts: 

* Estimation of the probabilities $p_ {ij}$ that gene *i* is mutated in sample *j* -assuming conditional independence between genes and samples-.
* Estimation of p-values using the Poisson-Binomial distribution, using the previous probabilities and the number of samples in which two genes are co-mutated. The corresponding null hypothesis $H_0$ is that the mutational status of both genes is independent of each other. 


The second step is the estimation of the p-values using these probabilities and the number of samples where two genes are co-altered. *Rediscover* offers different functions depending on the desired study:

* **`getMutex`** if the user wants to evaluate if genes are mutually exclusive.
* **`getMutexAB`** if the user wants to evaluate if genes are mutually exclusive with respect to another event (amplifications, deletions, etc...)
* **`getMutexGroup`** will be used when the user wants to obtain the probability that a certain group of genes being mutually exclusive. Unlike the `getMutex` function, in this case the users introduces the set of genes of interest. 


*Rediscover* also provides a function to integrate its usage with <a href="https://bioconductor.org/packages/release/bioc/html/maftools.html" target="_blank"> `maftools`</a> and <a href="https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html" target="_blank">`TCGAbiolinks`</a>. Specifically, we added to the function `somaticInteractions` from `maftools` our analyses based on the Poisson-Binomial distribution resulting in a new function called **`discoversomaticInteractions`**.
