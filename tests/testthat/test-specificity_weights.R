


load(system.file("testdata", "valid_data.RData", package = "LittleBites"))


test_that("spm output is numeric", {
  specificity <- max_spm(c('Alec' = 1, 'Alexander' = 2, 'Sarah' = 50, 'Harper' = 1, 'Sam' = 1, 'Jimmy' = 0))
  expect_true(is.numeric(specificity))
})

test_that("spm output is length 1", {
  specificity <- max_spm(c('Alec' = 1, 'Alexander' = 2, 'Sarah' = 50, 'Harper' = 1, 'Sam' = 1, 'Jimmy' = 0))
  expect_true(length(specificity)==1)
})

test_that("spm output from apply is a vector", {
  specificity <- apply(reference_valid, 1, max_spm)
  expect_true(is.vector(specificity))
})

test_that("spm output from apply is length of input", {
  specificity <- apply(reference_valid, 1, max_spm)
  expect_true(length(specificity) == nrow(reference_valid))
})


