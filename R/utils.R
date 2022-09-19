#' Creating a Class by Attribute Matrix
#'
#' Automating the creation of Class by Attribute Matrix
#'
#' @param natt An integer containing the number of assessed attributes.
#'
#' @return `profiles` A tibbler containing a class by attribute matrix listing
#' which attributes are mastered by each latent class.
#'
#' @export
#'
#' @examples
#' bin_profile(natt = 3)
bin_profile <- function(natt) {
  profiles <- rep(list(c(0L, 1L)), natt) %>%
    rlang::set_names(glue::glue("att_{seq_len(natt)}")) %>%
    expand.grid() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(total = rowSums(.)) %>%
    dplyr::select(tidyselect::everything(), .data$total) %>%
    dplyr::arrange(.data$total, -c(.data$total)) %>%
    dplyr::select(-.data$total)
  return(profiles)
}

#' Calculate the Number of Shards and Simultaneous Chains
#'
#' Calculating the number of shards and simultaneous chains.
#'
#' @param num_respondents An integer specifying the number of respondents.
#' @param num_responses An integer specifying the number of responses.
#' @param num_chains An integer specifying the number of chains that need to be
#' run.
#'
#' @return `ret` A list containing the number of shards to use within each chain
#' and the number of chains to run in parallel.
#'
#' @export
#'
#' @examples
#' shard_calculator(num_respondents = 1000, num_responses = 5000, num_chains = 4)
shard_calculator <- function(num_respondents, num_responses, num_chains) {
  max_shards <- parallel::detectCores() - 1

  possible_shards <- vector()
  possible_shards[1] <- 1
  kk <- 2

  for(jj in 2:max_shards) {
    if(num_respondents %% jj == 0 & num_responses %% jj == 0) {
      possible_shards[kk] <- jj
      kk <- kk + 1
    }
  }

  possible_parallel_chains <- floor(parallel::detectCores() / possible_shards)

  for(kk in 1:length(possible_parallel_chains)) {
    if(possible_parallel_chains[kk] >= parallel::detectCores()) {
      possible_parallel_chains[kk] <- parallel::detectCores() - 1
    }

    if(possible_parallel_chains[kk] > num_chains) {
      possible_parallel_chains[kk] <- num_chains
    }
  }

  optimal_config <- tibble::tibble(parallel_chains = possible_parallel_chains,
                                   threads_per_chain = possible_shards) %>%
    dplyr::mutate(total_cores = .data$parallel_chains * .data$threads_per_chain,
                  parallel_chains = dplyr::case_when(.data$total_cores >=
                                                       parallel::detectCores() ~
                                                       floor((parallel::detectCores() -
                                                                1) /
                                                               .data$threads_per_chain),
                                                     T ~ .data$parallel_chains),
                  total_cores = .data$parallel_chains *
                    .data$threads_per_chain) %>%
    dplyr::group_by(.data$parallel_chains) %>%
    dplyr::filter(.data$threads_per_chain == max(.data$threads_per_chain)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(num_sets = ceiling(num_chains / .data$parallel_chains)) %>%
    dplyr::group_by(.data$num_sets) %>%
    dplyr::filter(.data$threads_per_chain == max(.data$threads_per_chain)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$total_cores == max(.data$total_cores)) %>%
    dplyr::filter(.data$parallel_chains == max(.data$parallel_chains))

  ret <- list(n_shards_to_use = optimal_config$threads_per_chain[1],
              parallel_chains = optimal_config$parallel_chains[1])

  return(ret)
}
