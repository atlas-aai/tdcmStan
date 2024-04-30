test_that("test check_multi_attribute", {
  q_matrix <- tibble::tibble(att_1 = c(1, 0, 1, 0),
                             att_2 = c(0, 1, 0, 1))

  output <- tdcmStan:::check_multi_attribute(q_matrix)
  testthat::expect_equal(output, q_matrix)

  q_matrix <- tibble::tibble(att_1 = c(1, 0, 1, 0),
                             att_2 = c(0, 1, 1, 1))

  err <- rlang::catch_cnd(tdcmStan:::check_multi_attribute(q_matrix))
  expect_match(err$message, "Q-matrix cannot have items measuring multiple")
})
