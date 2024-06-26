test_that("Create TDCM works", {
  q_matrix <- tibble::tibble(att_1 = c(rep(1, 8), rep(0, 16)),
                             att_2 = c(rep(0, 8), rep(1, 8), rep(0, 8)),
                             att_3 = c(rep(0, 16), rep(1, 8)))

  tdcm_stan <- create_stan_tdcm(q_matrix)
  tdcm_stan %>%
    readr::write_lines(testthat::test_path("data", "tdcm-stan.stan"))

  tdcm_stan <- readr::read_lines(testthat::test_path("data", "tdcm-stan.stan"))
  true_tdcm_stan <-
    readr::read_lines(testthat::test_path("data", "madison-tdcm.stan"))

  testthat::expect_equal(tdcm_stan, true_tdcm_stan)
})

test_that("Create TDCM with multi-attribute items works", {
  q_matrix <- tibble::tibble(att_1 = c(rep(1, 4), rep(0, 4)),
                             att_2 = c(rep(0, 3), rep(1, 5)))

  tdcm_stan <- create_stan_tdcm(q_matrix)
  tdcm_stan %>%
    readr::write_lines(testthat::test_path("data", "multi-att-stan.stan"))

  tdcm_stan <- readr::read_lines(testthat::test_path("data",
                                                     "multi-att-stan.stan"))
  true_tdcm_stan <-
    readr::read_lines(testthat::test_path("data", "multi-att-tdcm.stan"))

  testthat::expect_equal(tdcm_stan, true_tdcm_stan)
})
