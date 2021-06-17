context('BOSO Wrapper')
library(BOSO)

x <- 

test_that('BOSO expects cplex API to be installed',{
  expect_equal(requireNamespace('cplexAPI', quietly = T), TRUE)
})

test_that('BOSO expects several inputs',{
  expect_error(BOSO())
  expect_error(BOSO(1))
  expect_error(BOSO(1,1))
  expect_error(BOSO(1,1))
  expect_error(BOSO(1,1,1))
})

test_that('BOSO expects several inputs',{
  expect_error(BOSO())
  expect_error(BOSO(1))
  expect_error(BOSO(1,1))
  expect_error(BOSO(1,1))
  expect_error(BOSO(1,1,1))
})