test_that("Create Fungible TDCM with no common items works", {
  q_matrix <- tibble::tibble(att_1 = rep(1, 24))

  fng_tdcm_no_common_stan <- create_fng_no_common_items_stan_tdcm(q_matrix)
  fng_tdcm_no_common_stan %>%
    readr::write_lines(testthat::test_path("data/fng-tdcm-no-common-stan.stan"))

  fng_tdcm_no_common_stan <-
    readr::read_lines(testthat::test_path("data/fng-tdcm-no-common-stan.stan"))
  true_fng_tdcm_no_common_stan <-
    readr::read_lines(testthat::test_path("data/dlm-fng-tdcm-no-common.stan"))

  testthat::expect_equal(fng_tdcm_no_common_stan, true_fng_tdcm_no_common_stan)
})
