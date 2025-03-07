


load(system.file("testdata", "valid_data.RData", package = "LittleBites"))
specificity <- apply(reference_valid |> log1p(), 1, LittleBites::max_spm)

### output should be a matrix
test_that("output is a dataframe", {
  result <- subtraction(bulk = bulk_valid,
                        reference = reference_valid,
                        cell_types_matrix = cell_types_matrix,
                        training_matrix = gt_train_valid,
                        sample_name_separator = 'r',
                        specificity_weights = specificity,
                        verbose = F)
  expect_true(is.data.frame(result))
})

test_that("output is numeric", {
  result <- subtraction(bulk = bulk_valid,
                        reference = reference_valid,
                        cell_types_matrix = cell_types_matrix,
                        training_matrix = gt_train_valid,
                        sample_name_separator = 'r',
                        specificity_weights = specificity,
                        verbose = F)
  expect_true(all(sapply(result, is.numeric)))
})

test_that("output has same dimensions as bulk input", {
  result <- subtraction(bulk = bulk_valid,
                        reference = reference_valid,
                        cell_types_matrix = cell_types_matrix,
                        training_matrix = gt_train_valid,
                        sample_name_separator = 'r',
                        specificity_weights = specificity,
                        verbose = F)
  expect_true(all(dim(bulk_valid) == dim(result)))
})

test_that("output of test data is non-negative", {
  result <- subtraction(bulk = bulk_valid,
                        reference = reference_valid,
                        cell_types_matrix = cell_types_matrix,
                        training_matrix = gt_train_valid,
                        sample_name_separator = 'r',
                        specificity_weights = specificity,
                        verbose = F)
  expect_true(all(result >= 0))
})


### subtraction is deterministic, so if randomness starts to come into play, something has gone wrong
test_that("function produces consistent results", {
  result1 <- subtraction(bulk = bulk_valid,
                         reference = reference_valid,
                         cell_types_matrix = cell_types_matrix,
                         training_matrix = gt_train_valid,
                         sample_name_separator = 'r',
                         specificity_weights = specificity,
                         verbose = Fs)
  result2 <- subtraction(bulk = bulk_valid,
                         reference = reference_valid,
                         cell_types_matrix = cell_types_matrix,
                         training_matrix = gt_train_valid,
                         sample_name_separator = 'r',
                         specificity_weights = specificity,
                         verbose = F)
  expect_equal(result1, result2)
})



bulk_valid_esep <- bulk_valid |>
  dplyr::rename_with(~ stringr::str_replace_all(.x, "r", "e"))
cell_types_matrix_esep <- cell_types_matrix
rownames(cell_types_matrix_esep) <- stringr::str_replace_all(rownames(cell_types_matrix_esep), "r", "e")

test_that("having a different separator in the cell_types_matrix and the bulk data should throw an error", {
  expect_error(subtraction(bulk = bulk_valid_esep,
                              reference = reference_valid,
                              cell_types_matrix = cell_types_matrix,
                              training_matrix = gt_train_valid,
                              specificity_weights = specificity,
                              sample_name_separator = 'e',
                              verbose = F))
})

test_that("so long as string separators match internally, the results should be the same", {
  result1 <- subtraction(bulk = bulk_valid,
                         reference = reference_valid,
                         cell_types_matrix = cell_types_matrix,
                         training_matrix = gt_train_valid,
                         sample_name_separator = 'r',
                         specificity_weights = specificity,
                         verbose = F)
  result2 <- subtraction(bulk = bulk_valid_esep,
                         reference = reference_valid,
                         cell_types_matrix = cell_types_matrix_esep,
                         training_matrix = gt_train_valid,
                         sample_name_separator = 'e',
                         specificity_weights = specificity,
                         verbose = F)
  expect_true(all(result1 == result2))
})

