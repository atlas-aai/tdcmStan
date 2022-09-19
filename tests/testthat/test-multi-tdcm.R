test_that("Create Multi-Threaded TDCM works", {
  q_matrix <- tibble::tibble(att_1 = c(rep(1, 8), rep(0, 16)),
                             att_2 = c(rep(0, 8), rep(1, 8), rep(0, 8)),
                             att_3 = c(rep(0, 16), rep(1, 8)))

  multi_tdcm_stan <- create_threaded_stan_tdcm(q_matrix)
  multi_tdcm_stan %>%
    readr::write_lines(testthat::test_path("data/multi-tdcm-stan.stan"))

  multi_tdcm_stan <-
    readr::read_lines(testthat::test_path("data/multi-tdcm-stan.stan"))
  true_multi_tdcm_stan <-
    readr::read_lines(testthat::test_path("data/automated-madison-tdcm.stan"))

  testthat::expect_equal(multi_tdcm_stan, true_multi_tdcm_stan)
})
