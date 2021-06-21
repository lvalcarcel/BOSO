context('BOSO Wrapper')
library(BOSO)

x <- matrix(rnorm(100*10), nrow = 100, ncol = 10)
y <- matrix(rnorm(100*1), nrow = 100, ncol = 1)
xval <- matrix(rnorm(100*10), nrow = 100, ncol = 10)
yval <- matrix(rnorm(100*1), nrow = 100, ncol = 1)

test_that('BOSO expects cplex API to be installed',{
  expect_equal(requireNamespace('cplexAPI', quietly = T), TRUE)
})

test_that('BOSO expects several inputs',{
  expect_error(BOSO())
  expect_error(BOSO(x))
  expect_error(BOSO(x,y))
  expect_error(BOSO(x,y,xval))
})

