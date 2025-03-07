

valid_expr <- c('gene1' = 0, 'gene2' = 10, 'gene3' = 3, 'gene4' = 8)
valid_gt <- c('gene1' = 0, 'gene2' = 1, 'gene3' = 0, 'gene4' = 1)

invalid_gt_na <- c('gene1' = 0, 'gene2' = 1, 'gene3' = NA, 'gene4' = 1)
invalid_gt_wrong_length <- c('gene1' = 0, 'gene2' = 1, 'gene4' = 1)

get_tpr(valid_expr, valid_gt, 4)

#load(system.file("testdata", "valid_data.RData", package = "LittleBites"))

## get_tpr

test_that("tpr output is numeric", {
  tpr <- get_tpr(valid_expr, valid_gt, 4)
  expect_true(is.numeric(tpr))
})

test_that("tpr output is 1 when threshold is 0", {
  tpr <- get_tpr(valid_expr, valid_gt, 0)
  expect_true(tpr == 1)
})

test_that("tpr output is 0 when threshold is infinite", {
  tpr <- get_tpr(valid_expr, valid_gt, Inf)
  expect_true(tpr == 0)
})


test_that("get_tpr error when given the wrong truth length", {
  expect_error(get_tpr(valid_expr, invalid_gt_wrong_length, 1), 'expression and truth vectors must be the same length')
})

### get_fpr

test_that("tpr output is numeric", {
  tpr <- get_tpr(valid_expr, valid_gt, 4)
  expect_true(is.numeric(tpr))
})

test_that("fpr output is 0 when threshold is infinite", {
  fpr <- get_fpr(valid_expr, valid_gt, Inf)
  expect_true(fpr == 0)
})

test_that("fpr output is 0 when threshold is infinite", {
  fpr <- get_fpr(valid_expr, valid_gt, Inf)
  expect_true(fpr == 0)
})

### get_fdr

test_that("fdr output is numeric when not infinite", {
  fdr <- get_fdr(valid_expr, valid_gt, 4)
  expect_true(is.numeric(fdr))
})

test_that("fdr output is NA when threshold is infinite", {
  fdr <- get_fdr(valid_expr, valid_gt, Inf)
  expect_true(is.na(fdr))
})

test_that("fdr output is NA when threshold is infinite", {
  fdr <- get_fdr(valid_expr, valid_gt, Inf)
  expect_true(is.na(fdr))
})



