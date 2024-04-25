library(tidyverse)
devtools::load_all()

qmatrix <- tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 1),
                          att_2 = c(0, 1, 0, 1, 1, 1))

create_stan_tdcm(q_matrix = qmatrix) %>%
  write_lines("inst/stan/tdcm.stan")

create_threaded_stan_tdcm(q_matrix = qmatrix) %>%
  write_lines("inst/stan/threaded-tdcm.stan")

create_fng_stan_tdcm(q_matrix = qmatrix) %>%
  write_lines("inst/stan/fng-tdcm.stan")

create_fng_no_common_item_tdcm(q_matrix = qmatrix) %>%
  write_lines("inst/stan/fng-no-common-tdcm.stan")
