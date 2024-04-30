test_that("bin_profile works", {
  output <- bin_profile(3)

  expected_qmatrix <- tibble::tibble(att_1 = c(0, 1, 0, 0, 1, 1, 0, 1),
                                     att_2 = c(0, 0, 1, 0, 1, 0, 1, 1),
                                     att_3 = c(0, 0, 0, 1, 0, 1, 1, 1))

  testthat::expect_equivalent(3, ncol(output))
  testthat::expect_equivalent(2^3, nrow(output))
  testthat::expect_equivalent(c("att_1", "att_2", "att_3"), colnames(output))
  testthat::expect_equivalent(output, expected_qmatrix)
})

test_that("shard_calculator works", {
  num_respondents <- 1000
  num_responses <- 5
  num_chains <- 4

  output <- shard_calculator(num_respondents, num_responses, num_chains)

  testthat::expect_equivalent(length(output), 2)
  testthat::expect_equivalent(length(output$n_shards_to_use), 1)
  testthat::expect_equivalent(length(output$parallel_chains), 1)
  testthat::expect_equivalent(output$n_shards_to_use %% 1, 0)
  testthat::expect_equivalent(output$parallel_chains %% 1, 0)
  testthat::expect_gt(output$n_shards_to_use, 0)
  testthat::expect_gt(output$parallel_chains, 0)
})
