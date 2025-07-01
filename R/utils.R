#' Creating a Class by Attribute Matrix
#'
#' Automating the creation of Class by Attribute Matrix
#'
#' @param natt An integer containing the number of assessed attributes.
#'
#' @return `profiles` A tibble containing a class by attribute matrix listing
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
    dplyr::select(tidyselect::everything(), "total") %>%
    dplyr::arrange(.data$total, -c(.data$total)) %>%
    dplyr::select(-"total")
  return(profiles)
}

#' Calculate the Number of Shards and Simultaneous Chains
#'
#' Calculating the number of shards and simultaneous chains.
#'
#' @param num_respondents An integer specifying the number of respondents.
#' @param num_responses An integer specifying the number of responses (i.e.,
#' the total number of items completed across all of the respondents).
#' @param num_chains An integer specifying the number of chains that need to be
#' run.
#'
#' @return `ret` A list containing the number of shards to use within each chain
#' and the number of chains to run in parallel.
#'
#' @export
#'
#' @examples
#' shard_calculator(num_respondents = 1000, num_responses = 5000,
#'                  num_chains = 4)
shard_calculator <- function(num_respondents, num_responses, num_chains) {
  max_shards <- parallel::detectCores() - 1

  possible_shards <- vector()
  possible_shards[1] <- 1
  kk <- 2

  for (jj in 2:max_shards) {
    if (num_respondents %% jj == 0 && num_responses %% jj == 0) {
      possible_shards[kk] <- jj
      kk <- kk + 1
    }
  }

  possible_parallel_chains <- floor(parallel::detectCores() / possible_shards)

  for (kk in seq_len(length(possible_parallel_chains))) {
    if (possible_parallel_chains[kk] >= parallel::detectCores()) {
      possible_parallel_chains[kk] <- parallel::detectCores() - 1
    }

    if (possible_parallel_chains[kk] > num_chains) {
      possible_parallel_chains[kk] <- num_chains
    }
  }

  optimal_config <- tibble::tibble(parallel_chains = possible_parallel_chains,
                                   threads_per_chain = possible_shards) %>%
    dplyr::mutate(total_cores = .data$parallel_chains * .data$threads_per_chain,
                  parallel_chains =
                    dplyr::case_when(.data$total_cores >=
                                       parallel::detectCores() ~
                                       floor((parallel::detectCores() - 1) /
                                               .data$threads_per_chain),
                                     TRUE ~ .data$parallel_chains),
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

calc_priors <- function(q_matrix, fungible) {
  if (fungible) {
    int0_priors <- glue::glue("l_0 ~ normal(0, 2);")
    mef_priors <- glue::glue("l_1 ~ lognormal(0, 1);")
    int2_priors <- ""
  } else {
    int0_priors <- glue::glue("l{1:nrow(q_matrix)}_0 ~ normal(0, 2);")

    mef_priors <- q_matrix %>%
      tibble::rowid_to_column("item_id") %>%
      tidyr::pivot_longer(cols = c(-"item_id"), names_to = "attr",
                          values_to = "meas") %>%
      dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                          "att_"))) %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::select(-"meas") %>%
      dplyr::mutate(param =
                      glue::glue("l{item_id}_1{attr} ~ lognormal(0, 1);")) %>%
      dplyr::pull(.data$param)

    aug_q_matrix <- q_matrix %>%
      dplyr::rowwise() %>%
      dplyr::mutate(total =
                      sum(dplyr::c_across(where(is.numeric)))) %>%
      tibble::rowid_to_column("item_id")

    multi_att_items <- aug_q_matrix %>%
      dplyr::filter(.data$total > 1)

    if (nrow(multi_att_items) == 0) {
      int2 <- ""
      int2_priors <- ""
    } else {
      int2 <- multi_att_items %>%
        dplyr::filter(.data$total == 2) %>%
        tidyr::pivot_longer(cols = c(-"item_id"), names_to = "attr",
                            values_to = "meas") %>%
        dplyr::filter(.data$meas == 1) %>%
        dplyr::group_by(.data$item_id) %>%
        dplyr::mutate(att_num = dplyr::row_number(),
                      att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                                 .data$att_num == 2 ~ "att2")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                            "att_"))) %>%
        dplyr::select(-"meas") %>%
        tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
        dplyr::mutate(param =
                        glue::glue("real<lower=-1 * fmin(l{item_id}_1{att1}, ",
                                   "l{item_id}_1{att2})> ",
                                   "l{item_id}_2{att1}{att2};")) %>%
        dplyr::pull(.data$param)
      int2_priors <- multi_att_items %>%
        dplyr::filter(.data$total == 2) %>%
        dplyr::select(-"total") %>%
        tidyr::pivot_longer(cols = c(-"item_id"), names_to = "attr",
                            values_to = "meas") %>%
        dplyr::filter(.data$meas == 1) %>%
        dplyr::group_by(.data$item_id) %>%
        dplyr::mutate(att_num = dplyr::row_number(),
                      att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                                 .data$att_num == 2 ~ "att2")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                            "att_"))) %>%
        dplyr::select(-"meas") %>%
        tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
        dplyr::mutate(param =
                        glue::glue("l{item_id}_2{att1}{att2} ~ ",
                                   "normal(0, 2);")) %>%
        dplyr::pull(.data$param)
    }
  }

  ret_list <- list(intercept = int0_priors,
                   main_effect = mef_priors,
                   interactions = int2_priors)

  return(ret_list)
}

calc_pi_mat <- function(q_matrix, profs, profile_item_interactions, fungible) {
  if (fungible) {
    pi_mat <- tibble::tibble(profile = rep(1:(2^ncol(q_matrix)),
                                           each = nrow(q_matrix)),
                             item_id = rep(seq_len(nrow(q_matrix)),
                                           times = (2^ncol(q_matrix)))) %>%
      dplyr::left_join(profs %>%
                         tibble::rowid_to_column("profile") %>%
                         tidyr::pivot_longer(cols = c(-"profile"),
                                             names_to = "att_mastered",
                                             values_to = "mastered"),
                       by = "profile", relationship = "many-to-many") %>%
      dplyr::left_join(q_matrix %>%
                         tibble::rowid_to_column("item_id") %>%
                         tidyr::pivot_longer(cols = c(-"item_id"),
                                             names_to = "att_measured",
                                             values_to = "measured"),
                       by = "item_id", relationship = "many-to-many") %>%
      dplyr::filter(.data$att_mastered == .data$att_measured) %>%
      dplyr::mutate(int0 = glue::glue("l_0"),
                    need_param = .data$mastered * .data$measured,
                    attribute = as.numeric(stringr::str_remove(.data$att_measured,
                                                               "att_")),
                    mef =
                      dplyr::case_when(.data$need_param == 0 ~ NA_character_,
                                       .data$need_param > 0 ~
                                         as.character(glue::glue("l_1"))
                      )) %>%
      dplyr::select(-"att_measured", -"attribute", -"measured",
                    -"mastered", -"need_param") %>%
      tidyr::pivot_wider(names_from = "att_mastered", values_from = "mef") %>%
      dplyr::left_join(profile_item_interactions %>%
                         dplyr::rename(int2 = "param"),
                       by = c("profile", "item_id"),
                       relationship = "many-to-many") %>%
      tidyr::unite(col = "param", c(-"profile", -"item_id"), sep = "+",
                   na.rm = TRUE) %>%
      dplyr::mutate(stan_pi =
                      as.character(glue::glue("pi[{item_id},{profile}] = ",
                                              "inv_logit({param});")))
  } else {
    pi_mat <- tibble::tibble(profile = rep(1:(2^ncol(q_matrix)),
                                           each = nrow(q_matrix)),
                             item_id = rep(seq_len(nrow(q_matrix)),
                                           times = (2^ncol(q_matrix)))) %>%
      dplyr::left_join(profs %>%
                         tibble::rowid_to_column("profile") %>%
                         tidyr::pivot_longer(cols = c(-"profile"),
                                             names_to = "att_mastered",
                                             values_to = "mastered"),
                       by = "profile", relationship = "many-to-many") %>%
      dplyr::left_join(q_matrix %>%
                         tibble::rowid_to_column("item_id") %>%
                         tidyr::pivot_longer(cols = c(-"item_id"),
                                             names_to = "att_measured",
                                             values_to = "measured"),
                       by = "item_id", relationship = "many-to-many") %>%
      dplyr::filter(.data$att_mastered == .data$att_measured) %>%
      dplyr::mutate(int0 = glue::glue("l{item_id}_0"),
                    need_param = .data$mastered * .data$measured,
                    attribute = as.numeric(stringr::str_remove(.data$att_measured,
                                                               "att_")),
                    mef =
                      dplyr::case_when(.data$need_param == 0 ~ NA_character_,
                                       .data$need_param > 0 ~
                                         as.character(glue::glue("l{item_id}_1",
                                                                 "{attribute}"))
                      )) %>%
      dplyr::select(-"att_measured", -"attribute", -"measured",
                    -"mastered", -"need_param") %>%
      tidyr::pivot_wider(names_from = "att_mastered", values_from = "mef") %>%
      dplyr::left_join(profile_item_interactions %>%
                         dplyr::rename(int2 = "param"),
                       by = c("profile", "item_id"),
                       relationship = "many-to-many") %>%
      tidyr::unite(col = "param", c(-"profile", -"item_id"), sep = "+",
                   na.rm = TRUE) %>%
      dplyr::mutate(stan_pi =
                      as.character(glue::glue("pi[{item_id},{profile}] = ",
                                              "inv_logit({param});")))
  }

  return(pi_mat)
}

# think it's done
parameter_block <- function(mod, t, int0, mef, int2) {
  if (mod == "TDCM") {
    if (all(int2 == "")) {
      stan_parameters <-
        glue::glue("parameters {{",
                   "  array[C, {t-1}] simplex[C] tau;",
                   "  simplex[C] Vc;",
                   glue::glue_collapse(glue::glue("  {int0}"), "\n"),
                   glue::glue_collapse(glue::glue("  {mef}"), "\n"),
                   "}}", .sep = "\n")
    } else {
      stan_parameters <-
        glue::glue("parameters {{",
                   "  array[C, {t-1}] simplex[C] tau;",
                   "  simplex[C] Vc;",
                   glue::glue_collapse(glue::glue("  {int0}"), "\n"),
                   glue::glue_collapse(glue::glue("  {mef}"), "\n"),
                   glue::glue_collapse(glue::glue("  {int2}"), "\n"),
                   "}}", .sep = "\n")
    }
  }

  return(stan_parameters)
}

# think it's done
gqs_block <- function(mod, t, multithread, repeated_items) {
  if (mod == "TDCM") {

    mod_lik <- calc_lik(mod, t, multithread, repeated_items)
    time_point_iters <- calc_time_loops(mod, t, mod_lik, stan_func = FALSE)
    time_point_iters <- stringr::str_replace(time_point_iters, "target \\+",
                                             "log_lik[j] ")
    time_point_iters <- stringr::str_replace(time_point_iters, "\\+\\=", "\\=")

    array_params <- stringr::str_flatten(rep("C", t), collapse = ", ")
    prob_transition_class_vector <- stringr::str_flatten(rep("C", t),
                                                         collapse = "*")
    prob_transition_class_vector <- stringr::str_c("J*",
                                                   prob_transition_class_vector)

    tau_declares <- tibble::tibble(time = 1:t) %>%
      dplyr::filter(.data$time != max(.data$time)) %>%
      dplyr::mutate(term =
                      dplyr::case_when(time == 1 ~
                                         glue::glue("  array[C] vector[J] tau_tmp;"),
                                       TRUE ~
                                         glue::glue("  array[C] vector[J] tau_tmp{time};")
                                       ))


    tau_declares <- glue::glue_collapse(tau_declares$term, sep = "\n")

    if (t == 2) {
      format_prob_transition_class_term <- ""
    } else {
      format_prob_transition_class_term <-
        stringr::str_c(",", stringr::str_flatten(rep('C', t-2), collapse = ','))
    }


    if (multithread) {
      intro <-  glue::glue("generated quantities {{",
                           "  vector[J] log_lik;",
                           "  array[J, {array_params}] real format_prob_transition_class;",
                           "  vector[{prob_transition_class_vector}] prob_transition_class;",
                           "  array[J] matrix[A, {t}] prob_resp_attr;",
                           # "  array[N, {t}] int<lower=0> y_rep;",
                           "  array[J, {t}] int j_class;",
                           "{tau_declares}",
                           .sep = "\n", .trim = FALSE)
    } else {
      intro <- glue::glue("generated quantities {{",
                          "  vector[J] log_lik;",
                          "  array[J, {array_params}] real prob_transition_class;",
                          "  array[J] matrix[A, {t}] prob_resp_attr;",
                          # "  array[N, {t}] int<lower=0> y_rep;",
                          "  array[J, {t}] int j_class;",
                          "{tau_declares}",
                          .sep = "\n", .trim = FALSE)
    }

    if (multithread) {
      lik <- "  log_lik = map_rect(person_loglik, beta, theta, xr, xi);"
    } else {
      lik <- glue::glue("  // Likelihood\n",
                        "  for (j in 1:J) {{\n",
                        "    matrix[C,C] tmp;\n",
                        "    array[{array_params}] real ps;\n",
                        "    {time_point_iters}",
                        .sep = "", .trim = FALSE)
    }

    mod_lik <- calc_lik(mod, t, multithread, repeated_items)
    time_point_iters <- calc_gqs_time_loops(mod, t, mod_lik)

    person_iter <- tibble::tibble(time = 1:t) %>%
      dplyr::mutate(c_term = glue::glue("c{time}"),
                    term = dplyr::case_when(.data$time != max(.data$time) ~
                                              glue::glue("(({c_term} - 1) * C^{t-time})"),
                                            TRUE ~ .data$c_term)) %>%
      dplyr::select("term")
    person_iter <- glue::glue_collapse(person_iter$term, sep = " + ")
    person_iter <- stringr::str_c(person_iter, glue::glue("((j - 1) * C^{t})"),
                                  sep = " + ")

    prob_params <- stringr::str_flatten(glue::glue("c{1:t}"), collapse = ",")

    gqs_buffer <- stringr::str_pad("", 4 + (2 * t), side = "left", pad = " ")

    contents <- glue::glue("<gqs_buffer>int iter = to_int(<person_iter>);\n",
                           "<gqs_buffer>format_prob_transition_class[j,<prob_params>] = prob_transition_class[iter];",
                           .trim = FALSE, .open = "<", .close = ">")

    gqs_loops <- time_loop(t, contents)

    if (multithread) {
      gqs_lik <- glue::glue("  prob_transition_class = map_rect(resp_transition, beta, theta, xr, xi);",
                            "",
                            "  for (j in 1:J) {{",
                            "<gqs_loops>",
                            "  }}",
                            .sep = "\n", .trim = FALSE, .open = "<",
                            .close = ">")
    } else {
      gqs_lik <- glue::glue("  // latent class probabilities",
                            "  for (j in 1:J) {",
                            "    array[<array_params>] real tmp;",
                            "    array [<array_params>] real prob_joint;",
                            "    <time_point_iters>",
                            "  }",
                            .sep = "\n", .open = "<", .close = ">",
                            .trim = FALSE)
    }

    # buffer <- stringr::str_pad("", 2 * t, side = "left", pad = " ")
    # buffer2 <- stringr::str_pad("", 2 * t + 2, side = "left", pad = " ")
    # prob_attr_vectors <-
    #   glue::glue_collapse(glue::glue("{buffer}vector[C] prob_attr_class_t{1:t};\n"),
    #                       sep = "\n")
    # prob_resp_attr <- tibble::tibble(time = 1:t) %>%
    #   dplyr::mutate(term = glue::glue("{buffer}prob_resp_attr[j,a,{time}] = sum(prob_attr_class_t{time});"))
    # prob_resp_attr <- glue::glue_collapse(prob_resp_attr$term, sep = "\n")
    #
    # prob_attr_class_t <- tibble::tibble(time = 1:t) %>%
    #   dplyr::mutate(filler = stringr::str_flatten(rep(",", t)),
    #                 filler = stringr::`str_sub<-`(.data$filler, .data$time,
    #                                               .data$time, value = "c"),
    #                 params = stringr::str_c("j,", .data$filler, sep = ""),
    #                 term = glue::glue("{buffer2}prob_attr_class_t{time}[c] = sum(to_matrix(format_prob_transition_class[{params}])) * Alpha[c,a];"))
    #
    # prob_attr_class_t <- glue::glue_collapse(prob_attr_class_t$term, sep = "\n")
    #
    # prob_resp_attr <- glue::glue("  for (j in 1:J) {{",
    #                              "    for (a in 1:A) {{",
    #                              "{prob_attr_vectors}",
    #                              "      for (c in 1:C) {{",
    #                              "{prob_attr_class_t}",
    #                              "      }}",
    #                              "{prob_resp_attr}",
    #                              "    }}",
    #                              "  }}",
    #                              .sep = "\n", .trim = FALSE)
    #
    # tau_tmp <- tibble::tibble(time = 1:t) %>%
    #   dplyr::mutate(tau_tmp_filler = stringr::str_flatten(rep(", ", t-1,
    #                                                           collapse = "")),
    #                 j_class = dplyr::case_when(.data$time == 1 ~
    #                                              glue::glue("    j_class[j, {time}] = categorical_rng(j_probs);"),
    #                                            .data$time == 2 ~
    #                                              glue::glue("    j_class[j, {time}] = categorical_rng(tau_tmp[j]);"),
    #                                            .data$time > 2 ~
    #                                              glue::glue("    j_class[j, {time}] = categorical_rng(tau_tmp{time-1}[j]);")),
    #                 tau_tmp = dplyr::case_when(.data$time == 1 ~
    #                                              glue::glue("    tau_tmp = tau[j_class[j, {time}]{tau_tmp_filler}];"),
    #                                            .data$time == max(.data$time) ~
    #                                              "",
    #                                            .data$time > 1 ~
    #                                              glue::glue("    tau_tmp{time} = tau[j_class[j, {time}]{tau_tmp_filler}];"))) %>%
    #   dplyr::select(-"time", -"tau_tmp_filler") %>%
    #   dplyr::mutate(term = stringr::str_c(.data$j_class, .data$tau_tmp,
    #                                       sep = "\n")) %>%
    #   dplyr::select("term")
    #
    # tau_tmp <- glue::glue_collapse(tau_tmp$term, sep = "\n")
    #
    # y_rep <- glue::glue("  for (j in 1:J) {",
    #                     "    vector[C] j_probs = Vc / sum(Vc);",
    #                     "<tau_tmp>",
    #                     "    for (t in 1:<t>) {",
    #                     "      for (m in 1:l[j, t]) {",
    #                     "        int i = ii[s[j, t] + m - 1, t];",
    #                     "        y_rep[s[j, t] + m - 1, t] = bernoulli_rng(pi[i, j_class[j, t]]);",
    #                     "      }",
    #                     "    }",
    #                     "  }",
    #                     .sep = "\n", .trim = FALSE, .open = "<",
    #                     .close = ">")
    #
    # gqs <- glue::glue("{prob_resp_attr}",
    #                   "",
    #                   "{y_rep}",
    #                   "}}", .sep = "\n", .trim = FALSE)

    # if (repeated_items) {
    #   gqs <- glue::glue("{prob_resp_attr",
    #                     "",
    #                     "  for (j in 1:J) {{",
    #                     "    vector[C] j_probs = Vc / sum(Vc);",
    #                     "    j_class[j, 1] = categorical_rng(j_probs);",
    #                     "    tau_tmp = tau[j_class[j, 1], ];",
    #                     "    j_class[j, 2] = categorical_rng(tau_tmp);",
    #                     "    for (t in 1:2) {{",
    #                     "      for (m in 1:l[j, t]) {{",
    #                     "        int i = ii[s[j, t] + m - 1, t];",
    #                     "        y_rep[s[j, t] + m - 1, t] = bernoulli_rng(pi[i, j_class[j, t]]);",
    #                     "      }}",
    #                     "    }}",
    #                     "  }}",
    #                     "}}", .sep = "\n", .trim = FALSE)
    # } else {
    #   gqs <- glue::glue("  for (j in 1:J) {{",
    #                     "    for (a in 1:A) {{",
    #                     "      vector[C] prob_attr_class_t1;",
    #                     "      vector[C] prob_attr_class_t2;",
    #                     "      for (c in 1:C) {{",
    #                     "        prob_attr_class_t1[c] = sum(prob_transition_class[j,c,]) * Alpha[c,a];",
    #                     "        prob_attr_class_t2[c] = sum(prob_transition_class[j,,c]) * Alpha[c,a];",
    #                     "      }}",
    #                     "      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);",
    #                     "      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);",
    #                     "    }}",
    #                     "  }}",
    #                     "",
    #                     "  for (j in 1:J) {{",
    #                     "    vector[C] j_probs = Vc / sum(Vc);",
    #                     "    j_class[j, 1] = categorical_rng(j_probs);",
    #                     "    tau_tmp = tau[j_class[j, 1], ];",
    #                     "    j_class[j, 2] = categorical_rng(tau_tmp);",
    #                     "    for (t in 1:2) {{",
    #                     "      for (m in 1:l[j, t]) {{",
    #                     "        int i = ii[s[j, t] + m - 1, t];",
    #                     "        y_rep[s[j, t] + m - 1, t] = bernoulli_rng(pi[i, j_class[j, t]]);",
    #                     "      }}",
    #                     "    }}",
    #                     "  }}",
    #                     "}}", .sep = "\n", .trim = FALSE)
    # }


    stan_generated_quantities <- glue::glue("{intro}",
                                            "",
                                            "{lik}",
                                            "",
                                            "{gqs_lik}",
                                            "}}",
                                            # "",
                                            # "{gqs}",
                                            .sep = "\n", .trim = FALSE)
  }

  return(stan_generated_quantities)
}

stan_func_block <- function(mod, t, multithread, repeated_items) {
  if (mod == "TDCM") {
    array_code <- tibble::tibble(time = 1:t) %>%
      dplyr::mutate(y_var = glue::glue("y{time}"),
                    ii_var = glue::glue("ii{time}"),
                    jj_var = glue::glue("jj{time}"),
                    s_var = glue::glue("s{time}"),
                    l_var = glue::glue("l{time}"),
                    array_iter = .data$time - 1,
                    y = glue::glue("    array[iis] int {y_var} = xi[((3 * {array_iter} * iis) + (2 * {array_iter} * ys) + 1):(((1 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys))];"),
                    ii = glue::glue("    array[iis] int {ii_var} = xi[(((1 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys) + 1):(((2 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys))];"),
                    jj = glue::glue("    array[iis] int {jj_var} = xi[(((2 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys) + 1):(((3 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys))];"),
                    s = glue::glue("    array[ys] int {s_var} = xi[(((3 + (3 * {array_iter})) * iis) + (2 * {array_iter} * ys) + 1):(((3 + (3 * {array_iter})) * iis) + ((1 + (2 * {array_iter})) * ys))];"),
                    l = glue::glue("    array[ys] int {l_var} = xi[(((3 + (3 * {array_iter})) * iis) + ((1 + (2 * {array_iter})) * ys) + 1):(((3 + (3 * {array_iter})) * iis) + ((2 + (2 * {array_iter})) * ys))];")) %>%
      tidyr::pivot_longer(cols = c("y", "ii", "jj", "s", "l"),
                          names_to = "var",
                          values_to = "stan_code") %>%
      dplyr::select("stan_code")

    ps_term <- glue::glue_collapse(
      stringr::str_sub(glue::glue("C{1:t}", sep = ""), 1, 1), sep = ", ")

    person_term <- glue::glue_collapse(
      stringr::str_sub(glue::glue("C{1:t}", sep = ""), 1, 1), sep = " * ")

    mod_lik <- calc_lik(mod, t, multithread, repeated_items)
    time_point_iters <- calc_time_loops(mod, t, mod_lik, stan_func = TRUE)

    lik <- glue::glue("    // Likelihood\n",
                      "    for (j in 1:J) {{\n",
                      "      matrix[C,C] tmp;\n",
                      "      {time_point_iters}",
                      .sep = "", .trim = FALSE)
    person_log_lik <- stringr::str_replace(lik, "person \\+\\=",
                                           "person[j] =")

    if (t == 2) {
      prob_trans_term <- ""
    } else {
      prob_trans_term <- rep(",C", t-2)
    }

    if (t == 2) {
      prob_joint_term <- ""
    } else {
      prob_joint_term <-
        glue::glue("array[{stringr::str_flatten(rep('C', t-2), collapse = ',')}] ",
                   .trim = FALSE)
    }

    resp_transition <- stringr::str_replace(lik, "// Likelihood",
                                                "// latent class probabilities")
    resp_transition <- stringr::str_replace(resp_transition,
                                            "      vector\\[C] tmp;\\n",
                                            glue::glue("      vector[C] tmp;\n      <prob_joint_term>matrix[C,C] prob_joint;\n",
                                                       .trim = FALSE, .open = "<",
                                                       .close = ">"))
    # resp_transition <- stringr::str_replace(resp_transition, "ps\\[",
    #                                         "prob_joint[")
    resp_transition <-
      stringr::str_replace(resp_transition,
                           "person \\+\\= log_sum_exp\\(tmp\\);",
                           "person[j] \\= ps;")

    log_sum_remove <- stringr::str_flatten(as.character(glue::glue('c{1:(t-1)}\\,')))

    resp_transition <- stringr::str_remove(
      resp_transition,
      glue::glue("      tmp\\[c1\\] \\+\\= log_sum_exp\\(ps\\[{log_sum_remove}\\]\\)\\;\\n"))

    # resp_transition <-
    #   stringr::str_replace(resp_transition,
    #                        "person \\+\\= log_sum_exp\\(tmp\\)\\;\\n",
    #                        "prob_transition_class[j] = exp(prob_joint) / sum(exp(prob_joint));\n")

    prob_class_term <- stringr::str_c("[j,",
                                      stringr::str_flatten(glue::glue("c{1:(t-2)}"),
                                            collapse = ","),
                                      stringr::str_flatten(rep(",", t-1),
                                                           collapse = ""),
                                      "]")
    classes_term <- stringr::str_flatten(glue::glue(",c{1:(t-2)}"),
                                         collapse = ",")

    resp_transition <-
      stringr::str_replace(resp_transition,
                           "prob_transition_class\\[j\\] \\= exp\\(prob_joint\\) / sum\\(exp\\(prob_joint\\)\\)\\;\\n    \\}",
                           glue::glue("prob_transition_class[j{classes_term}] = exp(prob_joint{prob_class_term}) / sum(exp(prob_joint{prob_class_term}));\n"))


    person_iter <- tibble::tibble(time = 1:t) %>%
      dplyr::mutate(c_term = glue::glue("c{time}"),
                    term = dplyr::case_when(.data$time != max(.data$time) ~
                                              glue::glue("(({c_term} - 1) * C^{t-time})"),
                                            TRUE ~ .data$c_term)) %>%
      dplyr::select("term")
    person_iter <- glue::glue_collapse(person_iter$term, sep = " + ")
    person_iter <- stringr::str_c(person_iter, glue::glue("((j - 1) * C^{t})"),
                                  sep = " + ")

    prob_params <- stringr::str_flatten(glue::glue("c{1:t}"), collapse = ",")

    multithread_buffer <- stringr::str_pad("", 4 + (2 * t), side = "left", pad = " ")

    contents <- glue::glue("{multithread_buffer}person[to_int({person_iter})] = prob_transition_class[j,{prob_params}];\n")

    multithread_text <- time_loop(t, contents)

    resp_transition <- stringr::str_remove(resp_transition,
                                           "person\\[j\\] \\= ps;")
    # resp_transition <- stringr::str_replace(resp_transition,
    #                                         "matrix\\[C,C\\] tmp;",
    #                                         glue::glue("matrix\\[C,C\\] tmp;\\\\n      vector[{person_term}] tmp_ps;"))

    person_vec_buffer <- stringr::str_flatten(rep("  ", t + 2))
    class_vec <- stringr::str_flatten_comma(glue::glue("c{1:t}"))
    person_iterator <- tibble::tibble(time = 1:t) %>%
      dplyr::arrange(desc(.data$time)) %>%
      dplyr::mutate(c = stringr::str_c("c", as.character(.data$time)),
                    c_iter = t - .data$time) %>%
      dplyr::mutate(c_term =
                      dplyr::case_when(.data$c_iter > 0 ~
                                         stringr::str_c(
                                           stringr::str_dup(" * C",
                                                            .data$c_iter)
                                         ),
                                       TRUE ~ ""),
                    c = dplyr::case_when(.data$time == max(.data$time) ~
                                           .data$c,
                                         TRUE ~ stringr::str_c("(",
                                                               .data$c,
                                                               " - 1)")),
                    term = stringr::str_c(.data$c, .data$c_term),
                    term = dplyr::case_when(.data$time == max(.data$time) ~
                                              .data$term,
                                            TRUE ~ stringr::str_c("(",
                                                                  .data$term,
                                                                  ")")))
    person_iterator <- stringr::str_flatten(person_iterator$term,
                                            collapse = " + ")
    person_vec_contents <- stringr::str_c(person_vec_buffer,
                                          "int cc = ",
                                          person_iterator,
                                          ";",
                                          "\n",
                                          person_vec_buffer,
                                          "tmp_ps[cc] = ps[",
                                          class_vec,
                                          "];",
                                          sep = "")
    person_vec_loop <- time_loop(t, person_vec_contents)
    person_vec_iterator <- glue::glue("(((j - 1) * {person_term}) + 1):(((j - 1) * {person_term}) + ({person_term}))")

    stan_functions <-
      glue::glue("functions {{\n",
                 "  real minmax (real x) {{\n",
                 "    if (x < .01) {{\n",
                 "      return 0.01;\n",
                 "    }}\n",
                 "\n",
                 "    if (x > 0.99) {{\n",
                 "      return 0.99;\n",
                 "    }}\n",
                 "\n",
                 "    return x;\n",
                 "  }}\n",
                 "\n",
                 "  vector sum_probs(vector beta, vector theta, array[] ",
                 "real xr, array[] int xi) {{\n",
                 "    int Z = num_elements(xi);\n",
                 "    int ys = xi[Z - 1];\n",
                 "    int iis = xi[Z];\n",
                 "\n",
                 glue::glue_collapse(array_code$stan_code, sep = "\n"),
                 "\n",
                 "    int I = xi[Z - 7];\n",
                 "    int N = xi[Z - 6];\n",
                 "    int C = xi[Z - 5];\n",
                 "    int A = xi[Z - 4];\n",
                 "    int J = xi[Z - 3];\n",
                 "    int M = xi[Z - 2];\n",
                 "\n",
                 "    vector[C] Vc = beta[1:C];\n",
                 "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];\n",
                 "    vector[C * C * @t-1?] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
                 "+ (C * C * @t-1?))];\n",
                 "\n",
                 "    array[@ps_term?] real ps;\n",
                 "    array[C, @t-1?, C] real tau_c;\n",
                 "    real person = 0;\n",
                 "\n",
                 "    matrix[I, C] pi_c;\n",
                 "    for(c in 1:C) {{\n",
                 "      for(i in 1:I) {{\n",
                 "        int ic = i + ((c-1) * I);\n",
                 "        pi_c[i, c] = pic[ic];\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "    for(tt in 1:@t-1?) {{\n",
                 "      for(c1 in 1:C) {{\n",
                 "        for(c2 in 1:C) {{\n",
                 "          int cc = c2 + ((c1 - 1) * C + (C * C * (tt - 1)));\n",
                 "          tau_c[c1, tt, c2] = tauc[cc];\n",
                 "        }}\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "@lik?",
                 "\n",
                 "    return [person]';\n",
                 "  }}\n",
                 "}}\n",
                 "\n",
                 "  vector person_loglik(vector beta, vector theta, array[] ",
                 "real xr, array[] int xi) {{\n",
                 "    int Z = num_elements(xi);\n",
                 "    int ys = xi[Z - 1];\n",
                 "    int iis = xi[Z];\n",
                 "\n",
                 glue::glue_collapse(array_code$stan_code, sep = "\n"),
                 "\n",
                 "    int I = xi[Z - 7];\n",
                 "    int N = xi[Z - 6];\n",
                 "    int C = xi[Z - 5];\n",
                 "    int A = xi[Z - 4];\n",
                 "    int J = xi[Z - 3];\n",
                 "    int M = xi[Z - 2];\n",
                 "\n",
                 "    vector[C] Vc = beta[1:C];\n",
                 "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];\n",
                 "    vector[C * C * @t-1?] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
                 "+ (C * C * @t-1?))];\n",
                 "\n",
                 "    array[@ps_term?] real ps;\n",
                 "    array[C, @t-1?, C] real tau_c;\n",
                 "    vector[J] person;\n",
                 "\n",
                 "    matrix[I, C] pi_c;\n",
                 "    for(c in 1:C) {{\n",
                 "      for(i in 1:I) {{\n",
                 "        int ic = i + ((c-1) * I);\n",
                 "        pi_c[i, c] = pic[ic];\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "    for(tt in 1:@t-1?) {{\n",
                 "      for(c1 in 1:C) {{\n",
                 "        for(c2 in 1:C) {{\n",
                 "          int cc = c2 + ((c1 - 1) * C + (C * C * (tt - 1)));\n",
                 "          tau_c[c1, tt, c2] = tauc[cc];\n",
                 "        }}\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "@person_log_lik?",
                 "\n",
                 "    return person;\n",
                 "  }}\n",
                 "}}\n",
                 "\n",
                 "  vector resp_transition(vector beta, vector theta, array[] ",
                 "real xr, array[] int xi) {{\n",
                 "    int Z = num_elements(xi);\n",
                 "    int ys = xi[Z - 1];\n",
                 "    int iis = xi[Z];\n",
                 "\n",
                 glue::glue_collapse(array_code$stan_code, sep = "\n"),
                 "\n",
                 "    int I = xi[Z - 7];\n",
                 "    int N = xi[Z - 6];\n",
                 "    int C = xi[Z - 5];\n",
                 "    int A = xi[Z - 4];\n",
                 "    int J = xi[Z - 3];\n",
                 "    int M = xi[Z - 2];\n",
                 "\n",
                 "    vector[C] Vc = beta[1:C];\n",
                 "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];\n",
                 "    vector[C * C * @t-1?] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
                 "+ (C * C * @t-1?))];\n",
                 "\n",
                 "    array[@ps_term?] real ps;\n",
                 "    array[C, @t-1?, C] real tau_c;\n",
                 "    array[J@prob_trans_term?] matrix[C,C] prob_transition_class;\n",
                 "\n",
                 "    vector[J * @person_term?] person;\n",
                 "\n",
                 "    matrix[I, C] pi_c;\n",
                 "    for(c in 1:C) {{\n",
                 "      for(i in 1:I) {{\n",
                 "        int ic = i + ((c-1) * I);\n",
                 "        pi_c[i, c] = pic[ic];\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "    for(tt in 1:@t-1?) {{\n",
                 "      for(c1 in 1:C) {{\n",
                 "        for(c2 in 1:C) {{\n",
                 "          int cc = c2 + ((c1 - 1) * C + (C * C * (tt - 1)));\n",
                 "          tau_c[c1, tt, c2] = tauc[cc];\n",
                 "        }}\n",
                 "      }}\n",
                 "    }}\n",
                 "\n",
                 "@resp_transition?",
                 "    vector[@person_term?] tmp_ps;\n",
                 "    @person_vec_loop?",
                 # "\n",
                 # "@multithread_text?",
                 "\n",
                 "    person[@person_vec_iterator?] = tmp_ps;\n",
                 "    }}\n",
                 "    return person;\n",
                 "  }}\n",
                 "}}\n", .sep = "", .open = "@", .close = "?")
  }

  return(stan_functions)
}

# think it's done
trans_data_block <- function(mod, t) {
  xi_code <- tibble::tibble(time = 1:t) %>%
    dplyr::mutate(array_iter = .data$time - 1,
                  y = glue::glue("    xi[i, ((3 * iis * {array_iter}) + (2 * ys * {array_iter}) + 1):(((1 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}))] = y[iilower:iiupper, {time}];"),
                  ii = glue::glue("    xi[i, (((1 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}) + 1):(((2 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}))] = ii[iilower:iiupper, {time}];"),
                  jj = glue::glue("    xi[i, (((2 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}) + 1):(((3 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}))] = jj[iilower:iiupper, {time}];"),
                  s = glue::glue("    xi[i, (((3 + (3 * {array_iter})) * iis) + (2 * ys * {array_iter}) + 1):(((3 + (3 * {array_iter})) * iis) + ((1 + (2 * {array_iter})) * ys))] = s[1:ys, {time}];"),
                  l = glue::glue("    xi[i, (((3 + (3 * {array_iter})) * iis) + ((1 + (2 * {array_iter})) * ys) + 1):(((3 + (3 * {array_iter})) * iis) + ((2 + (2 * {array_iter})) * ys))] = l[ylower:yupper, {time}];"),
                  I = dplyr::case_when(.data$time == max(.data$time) ~
                                         glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 1)] = I;"),
                                       TRUE ~ ""),
                  N = dplyr::case_when(.data$time == max(.data$time) ~
                                         glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 2)] = N / n_shards;"),
                                       TRUE ~ ""),
                  C = dplyr::case_when(.data$time == max(.data$time) ~
                                         glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 3)] = C;"),
                                       TRUE ~ ""),
                  A = dplyr::case_when(.data$time == max(.data$time) ~
                                         glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 4)] = A;"),
                                       TRUE ~ ""),
                  J = dplyr::case_when(.data$time == max(.data$time) ~
                                         glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 5)] = J / n_shards;"),
                                       TRUE ~ ""),
                  iis = dplyr::case_when(.data$time == max(.data$time) ~
                                           glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 6)] = iis;"),
                                         TRUE ~ ""),
                  ys = dplyr::case_when(.data$time == max(.data$time) ~
                                          glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 7)] = ys;"),
                                        TRUE ~ ""),
                  iis2 = dplyr::case_when(.data$time == max(.data$time) ~
                                            glue::glue("    xi[i, ((3 * {time} * iis) + (2 * {time} * ys) + 8)] = iis;"),
                                          TRUE ~ "")) %>%
    tidyr::pivot_longer(cols = -c("time", "array_iter"), names_to = "var",
                                  values_to = "stan_code") %>%
    dplyr::filter(.data$stan_code != "") %>%
    dplyr::select("stan_code")

  stan_transformed_data <-
    glue::glue("transformed data {{\n",
               "  int ys = num_elements(s) / <t> / n_shards;\n",
               "  int iis = num_elements(ii) / <t> / n_shards;\n",
               "\n",
               "  int M = iis;\n",
               "\n",
               "  array[n_shards, (2 * <t> * ys) + (3 * <t> * iis) + 8] int xi;\n",
               "\n",
               "  // an empty set of per-shard parameters\n",
               "  array[n_shards] vector[0] theta;\n",
               "\n",
               "  array[n_shards,1] real xr;\n",
               "  for(kk in 1:n_shards) {{\n",
               "    xr[kk, 1] = 1.0;\n",
               "  }}\n",
               "\n",
               "  // split into shards\n",
               "  for (i in 1:n_shards) {{\n",
               "    int ylower;\n",
               "    int yupper;\n",
               "    int iilower;\n",
               "    int iiupper;\n",
               "\n",
               "    ylower = ((i - 1) * ys) + 1;\n",
               "    yupper = i * ys;\n",
               "    iilower = ((i - 1) * iis) + 1;\n",
               "    iiupper = i * iis;\n",
               "\n",
               glue::glue_collapse(xi_code$stan_code, sep = "\n"),
               "\n",
               "  }}\n",
               "}}\n", .sep = "", .open = "<", .close = ">")

  return(stan_transformed_data)
}

# think it's done
trans_param_block <- function(mod, pi_mat, multithread) {
  if (mod == "TDCM") {
    if (multithread) {
      stan_transformed_parameters <-
        glue::glue("transformed parameters {{",
                   "  matrix[I,C] pi;",
                   "",
                   glue::glue_collapse(glue::glue("  {pi_mat$stan_pi}"), "\n"),
                   "",
                   "  array[I * C] real pic;",
                   "  for(c in 1:C) {{",
                   "    for(i in 1:I) {{",
                   "      int ic = i + ((c - 1) * I);",
                   "      pic[ic] = pi[i, c];",
                   "    }}",
                   "  }}",
                   "",
                   "  array[C * C * <t-1>] real tauc;",
                   "  for(tt in 1:<t-1>) {{",
                   "    for(c1 in 1:C) {{",
                   "      for(c2 in 1:C) {{",
                   "        int cc = c2 + ((c1 - 1) * C + (C * C * (tt - 1)));",
                   "        tauc[cc] = tau[c1, tt, c2];",
                   "      }}",
                   "    }}",
                   "  }}",
                   "",
                   "  // a set of shared parameters",
                   "  vector[C + (I * C) + (C * C * <t-1>)] beta;",
                   "  beta[1:C] = Vc[1:C];",
                   "  beta[(C + 1):(C + (I * C))] = to_vector(pic[1:(I * C)]);",
                   "  beta[(C + (I * C) + 1):(C + (I * C) + (C * C * <t-1>))] = to_vector(tauc[1:(C * C * <t-1>)]);",
                   "}}", .sep = "\n", .open = "<", .close = ">")
    } else {
      stan_transformed_parameters <-
        glue::glue("transformed parameters {{",
                   "  matrix[I,C] pi;",
                   "",
                   glue::glue_collapse(glue::glue("  {pi_mat$stan_pi}"), "\n"),
                   "}}", .sep = "\n")
    }
  }

  return(stan_transformed_parameters)
}

# think it's done
data_block <- function(mod, t, multithread) {
  if (mod == "TDCM") {
    if (multithread) {
      stan_data <-
        glue::glue("data {{",
                   "  int<lower=1> I;",
                   "  int<lower=1> J;",
                   "  int<lower=1> N;",
                   "  int<lower=1> C;",
                   "  int<lower=1> A;",
                   "  array[N, {t}] int<lower=1,upper=I> ii;",
                   "  array[N, {t}] int<lower=1,upper=J> jj;",
                   "  array[N, {t}] int<lower=0,upper=1> y;",
                   "  array[J, {t}] int<lower=1,upper=N> s;",
                   "  array[J, {t}] int<lower=1,upper=I> l;",
                   "  matrix[C,A] Alpha;",
                   "  int<lower=1> n_shards;",
                   "}}", .sep = "\n")
    } else {
      stan_data <-
        glue::glue("data {{",
                   "  int<lower=1> I;",
                   "  int<lower=1> J;",
                   "  int<lower=1> N;",
                   "  int<lower=1> C;",
                   "  int<lower=1> A;",
                   "  array[N, {t}] int<lower=1,upper=I> ii;",
                   "  array[N, {t}] int<lower=0,upper=1> y;",
                   "  array[J, {t}] int<lower=1,upper=N> s;",
                   "  array[J, {t}] int<lower=1,upper=I> l;",
                   "  matrix[C,A] Alpha;",
                   "}}", .sep = "\n")
    }
  }

  return(stan_data)
}

# think it's done
calc_lik <- function(mod, t, multithread, repeated_items) {
  if (multithread) {
    terms <- tibble::tibble(t = 1:t) %>%
      dplyr::mutate(term = glue::glue("y{t}[s{t}[j] + m - 1] * log(pi_c[i,c{t}]) + (1 - y{t}[s{t}[j] + m - 1]) * log(1 - pi_c[i,c{t}])"))

    lik <- stringr::str_c("log_items[m] = ",
                          glue::glue_collapse(terms$term, sep = " + "),
                          sep = "")
  } else if (!repeated_items) {
    buffer <- stringr::str_c("    ",
                             stringr::str_pad("", t * 2, side = "right",
                                              pad = " "))

    terms <- tibble::tibble(t = 1:t) %>%
      dplyr::mutate(term = glue::glue("<buffer>if(y[s[j, <t>] + m - 1, <t>] != 9) {{tmp<t> = y[s[j, <t>] + m - 1, <t>] * log(pi[i,c<t>]) + (1 - y[s[j, <t>] + m - 1, <t>]) * log(1 - pi[i,c<t>]);}}",
                                      .open = "<", .close = ">"))

    log_items <- glue::glue_collapse(glue::glue("tmp{1:t}"), sep = " + ")
    lik <- glue::glue("{buffer}real tmp{1:t} = 0;",
                      glue::glue_collapse(terms$term, sep = "\n"),
                      "{buffer}log_items[m] = {log_items};",
                      sep = "\n")
  } else {
    buffer <- stringr::str_c("    ",
                             stringr::str_pad("", t * 2, side = "right",
                                              pad = " "))

    terms <- tibble::tibble(t = 1:t) %>%
      dplyr::mutate(term = glue::glue("y[s[j, {t}] + m - 1, {t}] * log(pi[i,c{t}]) + (1 - y[s[j, {t}] + m - 1, {t}]) * log(1 - pi[i,c{t}])"))

    lik <- stringr::str_c(glue::glue("{buffer}log_items[m] = "),
                          glue::glue_collapse(terms$term, sep = " + "),
                          sep = "")
  }

  return(lik)
}

# think it's done
calc_time_loops <- function(mod, t, mod_lik, stan_func) {
  if (mod == "TDCM") {
    loops <- tibble::tibble(time = 1:t) %>%
      dplyr::mutate(buffer = stringr::str_pad("", width = 2 * .data$time + 2,
                                              side = "right", pad = " "),
                    term = glue::glue("<buffer>for (c<time> in 1:C) {\n",
                                      .open = "<", .close = ">"),
                    term = as.character(.data$term))

    buffer <- stringr::str_pad("", width = 2 * t + 4, side = "right", pad = " ")
    buffer5 <- stringr::str_pad("", width = 2 * t + 6, side = "right",
                                pad = " ")
    buffer6 <- stringr::str_pad("", width = 2 * t + 2, side = "right",
                                pad = " ")
    buffer2 <- stringr::str_pad("", width = 2 * t, side = "right", pad = " ")
    buffer3 <- stringr::str_pad("", width = 2 * t - 2, side = "right",
                                pad = " ")
    buffer4 <- stringr::str_pad("", width = 2 * t - 4, side = "right",
                                pad = " ")

    taus <- tidyr::crossing(combo1 = 1:t,
                            combo2 = 1:t) %>%
      dplyr::filter(.data$combo1 < .data$combo2) %>%
      dplyr::filter(.data$combo1 + 1 == .data$combo2) %>%
      dplyr::mutate(term =
                      glue::glue("log(tau[c{combo1}, {combo1}, c{combo2}])"))

    taus <- glue::glue_collapse(taus$term, sep = " + ")

    if (stan_func) {
      taus <- stringr::str_replace_all(taus, "tau", "tau_c")
    }

    ps_terms <- glue::glue_collapse(glue::glue("c{1:t}", sep = ""), sep = ", ")

    ps_sum_term <- stringr::str_flatten(as.character(glue::glue("c{1:(t-1)},")))

    time_point_iters <- glue::glue(glue::glue_collapse(loops$term, sep = "\n"),
                                   "\n",
                                   ifelse(stan_func,
                                          "<buffer>array[l1[j]] real log_items;\n",
                                          "<buffer>array[l[j, 1]] real log_items;\n"),
                                   ifelse(stan_func,
                                          "<buffer>for (m in 1:l1[j]) {\n",
                                          "<buffer>for (m in 1:l[j, 1]) {\n"),
                                   ifelse(stan_func,
                                          "<buffer5>int i = ii1[s1[j] + m - 1];\n",
                                          "<buffer5>int i = ii[s[j, 1] + m - 1, 1];\n"),
                                   "  <mod_lik>;\n",
                                   "<buffer>}\n",
                                   "<buffer>ps[<ps_terms>] = log(Vc[c1]) + <taus> + sum(log_items);\n",
                                   "<buffer6>}\n",
                                   "<buffer6>tmp[c1,c2] = log_sum_exp(ps[<ps_sum_term>]);\n",
                                   "<buffer2>}\n",
                                   "<buffer3>}\n",
                                   ifelse(stan_func,
                                          "<buffer3>person += log_sum_exp(tmp);\n",
                                          "<buffer3>target += log_sum_exp(tmp);\n"),
                                   ifelse(stan_func, "", "<buffer4>}\n"),
                                   .sep = "", trim = FALSE, .open = "<",
                                   .close = ">")

    return(time_point_iters)
  }
}

# think it's done
calc_gqs_time_loops <- function(mod, t, mod_lik) {
  if (mod == "TDCM") {
    loops <- tibble::tibble(time = 1:t) %>%
      dplyr::mutate(buffer = stringr::str_pad("", width = 2 * .data$time + 2,
                                              side = "right", pad = " "),
                    term = glue::glue("<buffer>for (c<time> in 1:C) {\n",
                                      .open = "<", .close = ">"),
                    term = as.character(.data$term))

    buffer <- stringr::str_pad("", width = 2 * t + 4, side = "right", pad = " ")
    buffer5 <- stringr::str_pad("", width = 2 * t + 6, side = "right",
                                pad = " ")
    buffer6 <- stringr::str_pad("", width = 2 * t + 2, side = "right",
                                pad = " ")
    buffer2 <- stringr::str_pad("", width = 2 * t, side = "right", pad = " ")
    buffer3 <- stringr::str_pad("", width = 2 * t - 2, side = "right",
                                pad = " ")
    buffer4 <- stringr::str_pad("", width = 2 * t - 4, side = "right",
                                pad = " ")

    taus <- tidyr::crossing(combo1 = 1:t,
                            combo2 = 1:t) %>%
      dplyr::filter(.data$combo1 < .data$combo2) %>%
      dplyr::filter(.data$combo1 + 1 == .data$combo2) %>%
      dplyr::mutate(term =
                      glue::glue("log(tau[c{combo1}, {combo1}, c{combo2}])"))

    taus <- glue::glue_collapse(taus$term, sep = " + ")

    ps_terms <- glue::glue_collapse(glue::glue("c{1:t}", sep = ""), sep = ", ")

    ps_sum_term <- stringr::str_c("c1",
                                  stringr::str_flatten(rep(",", t-1)), sep = "")

    time_point_iters <- glue::glue(glue::glue_collapse(loops$term, sep = "\n"),
                                   "\n",
                                   "<buffer>array[l[j, 1]] real log_items;\n",
                                   "<buffer>for (m in 1:l[j, 1]) {\n",
                                   "<buffer5>int i = ii[s[j, 1] + m - 1, 1];\n",
                                   "  <mod_lik>;\n",
                                   "<buffer>}\n",
                                   "<buffer>prob_joint[<ps_terms>] = log(Vc[c1]) + <taus> + sum(log_items);\n",
                                   "<buffer6>}\n",
                                   "<buffer2>}\n",
                                   "<buffer3>}\n",
                                   "<buffer3>prob_transition_class[j] = prob_joint;\n",
                                   .sep = "", trim = FALSE, .open = "<",
                                   .close = ">")

    return(time_point_iters)
  }
}

model_block <- function(mod, t, priors, multithread, repeated_items) {
  int0_priors <- priors$intercept
  mef_priors <- priors$main_effect
  int2_priors <- priors$interactions

  mod_lik <- calc_lik(mod, t, multithread, repeated_items)
  time_point_iters <- calc_time_loops(mod, t, mod_lik, stan_func = FALSE)

  if (mod == "TDCM") {
    ps_term <- glue::glue_collapse(
      stringr::str_sub(glue::glue("C{1:t}", sep = ""), 1, 1), sep = ", ")
    intro <-  glue::glue("model {{",
                         "  array[<ps_term>] real ps;",
                         .sep = "\n", .trim = FALSE, .open = "<", .close = ">")
    if (all(int2_priors == "")) {
      priors <- glue::glue("  // Priors",
                           glue::glue_collapse(glue::glue("  {int0_priors}"), "\n"),
                           glue::glue_collapse(glue::glue("  {mef_priors}"), "\n"),
                           .sep = "\n", .trim = FALSE)
    } else {
      priors <- glue::glue("  // Priors",
                           glue::glue_collapse(glue::glue("  {int0_priors}"), "\n"),
                           glue::glue_collapse(glue::glue("  {mef_priors}"), "\n"),
                           glue::glue_collapse(glue::glue("  {int2_priors}"), "\n"),
                           .sep = "\n", .trim = FALSE)
    }

    if (multithread) {
      lik <- "  target += sum(map_rect(sum_probs, beta, theta, xr, xi));"
    } else {
      lik <- glue::glue("  // Likelihood\n",
                        "  for (j in 1:J) {{\n",
                        "    matrix[C,C] tmp;\n",
                        "    {time_point_iters}",
                        .sep = "", .trim = FALSE)
    }

    stan_model <- glue::glue("{intro}",
                             "",
                             "{priors}",
                             "",
                             "{lik}",
                             "}}", .sep = "\n", .trim = FALSE)
  }

  return(stan_model)
}

time_loop <- function(t, contents) {
  time_loops_text <- tibble::tibble(time = 1:t) %>%
    dplyr::mutate(buffer = stringr::str_pad("", width = 2 * .data$time + 2,
                                            side = "right", pad = " "),
                  term = glue::glue("<buffer>for (c<time> in 1:C) {{\n",
                                    .open = "<", .close = ">"),
                  term = as.character(.data$term))
  time_loops_text <- glue::glue_collapse(time_loops_text$term, sep = "\n")

  closing_parentheses <- tibble::tibble(time = 1:t) %>%
    dplyr::mutate(buffer = stringr::str_pad("", width = 2 * .data$time + 2,
                                            side = "right", pad = " "),
                  term = glue::glue("<buffer>}}\n",
                                    .open = "<", .close = ">"),
                  term = as.character(.data$term)) %>%
    dplyr::arrange(dplyr::desc(.data$time))
  closing_parentheses <- glue::glue_collapse(closing_parentheses$term,
                                             sep = "\n")

  return_text <- glue::glue("{time_loops_text}",
                            "{contents}",
                            "{closing_parentheses}",
                            .sep = "\n")

  return(return_text)
}
