test_that("Create Fungible TDCM works", {
  q_matrix <- tibble::tibble(att_1 = rep(1, 24))

  fng_tdcm_stan <- create_fng_stan_tdcm(q_matrix)
  fng_tdcm_stan %>%
    readr::write_lines(testthat::test_path("data/fng-tdcm-stan.stan"))

  fng_tdcm_stan <-
    readr::read_lines(testthat::test_path("data/fng-tdcm-stan.stan"))
  true_fng_tdcm_stan <-
    readr::read_lines(testthat::test_path("data/dlm-fng-tdcm.stan"))

  testthat::expect_equal(fng_tdcm_stan, true_fng_tdcm_stan)
})

# test_that("Create Fungible TDCM with multi-attribute items works", {
#   q_matrix <- tibble::tibble(att_1 = c(rep(1, 4), rep(0, 4)),
#                              att_2 = c(rep(0, 3), rep(1, 5)))
#
#   tdcm_stan <- create_fng_stan_tdcm(q_matrix)
#   tdcm_stan %>%
#     readr::write_lines(testthat::test_path("data", "fng-multi-att-stan.stan"))
#
#   tdcm_stan <- readr::read_lines(
#     testthat::test_path("data",
#                         "fng-multi-att-stan.stan"))
#   true_tdcm_stan <-
#     readr::read_lines(testthat::test_path("data", "fng-multi-att-tdcm.stan"))
#
#   testthat::expect_equal(tdcm_stan, true_tdcm_stan)
# })
