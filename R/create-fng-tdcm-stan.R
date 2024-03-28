#' Creating Fungible TDCM Stan Code
#'
#' Automating the creation of fungible Stan code for a TDCM.
#'
#' @param q_matrix A tibble containing the assessment Q-matrix.
#'
#' @return `stan_code` A list containing the text for the Stan code blocks.
#'
#' @export
#'
#' @examples
#' qmatrix = tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 1), att_2 = c(0, 1, 0, 1, 1, 1))
#' create_fng_stan_tdcm(q_matrix = qmatrix)
create_fng_stan_tdcm <- function(q_matrix) {
  profs <- bin_profile(ncol(q_matrix))

  colnames(q_matrix) <- glue::glue("att_{1:ncol(q_matrix)}")

  int0 <- glue::glue("real l_0;")
  int0_priors <- glue::glue("l_0 ~ normal(0, 2);")

  mef <- glue::glue("real<lower=0> l_1;")
  mef_priors <- glue::glue("l_1 ~ lognormal(0, 1);")

  aug_q_matrix <- q_matrix %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(dplyr::c_across(where(is.numeric)))) %>%
    tibble::rowid_to_column("item_id")

  multi_att_items <- aug_q_matrix %>%
    dplyr::filter(.data$total > 1)

  if (nrow(multi_att_items) == 0) {
    int2 <- ""
    int2_priors <- ""
  } else {
    int2 <- multi_att_items %>%
      dplyr::filter(.data$total == 2) %>%
      tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "attr",
                          values_to = "meas") %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_num = dplyr::row_number(),
                    att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                               .data$att_num == 2 ~ "att2")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                          "att_"))) %>%
      dplyr::select(-.data$meas) %>%
      tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
      dplyr::mutate(param = glue::glue("real<lower=-1 * fmin(l{item_id}_1{att1}, l{item_id}_1{att2})> l{item_id}_2{att1}{att2};")) %>%
      dplyr::pull(.data$param)
    int2_priors <- multi_att_items %>%
      dplyr::filter(.data$total == 2) %>%
      dplyr::select(-.data$total) %>%
      tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "attr",
                          values_to = "meas") %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_num = dplyr::row_number(),
                    att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                               .data$att_num == 2 ~ "att2")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                          "att_"))) %>%
      dplyr::select(-.data$meas) %>%
      tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
      dplyr::mutate(param =
                      glue::glue("l{item_id}_2{att1}{att2} ~ normal(0, 2);")) %>%
      dplyr::pull(.data$param)
  }

  multi_item_q_matrix <- q_matrix %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(dplyr::c_across(where(is.numeric)))) %>%
    tibble::rowid_to_column("item_id") %>%
    dplyr::filter(.data$total == 2) %>%
    dplyr::select(-.data$total)

  if (nrow(multi_item_q_matrix) > 0) {
    items_with_interactions <- multi_item_q_matrix %>%
      tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "att",
                          values_to = "meas") %>%
      dplyr::mutate(meas_att = as.numeric(stringr::str_remove(.data$att,
                                                              "att_"))) %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_row = dplyr::row_number(),
                    att_row = stringr::str_c("att_",
                                             as.character(.data$att_row))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$att, -.data$meas) %>%
      tidyr::pivot_wider(names_from = "att_row", values_from = "meas_att") %>%
      dplyr::mutate(param =
                      as.character(glue::glue("l{item_id}_2{att_1}{att_2}"))) %>%
      dplyr::select(.data$item_id, .data$param)

    profile_item_interactions <-
      tibble::tibble(profile = rep(1:(2^ncol(q_matrix)),
                                   each = nrow(q_matrix)),
                     item_id = rep(seq_len(nrow(q_matrix)),
                                   times = (2^ncol(q_matrix)))) %>%
      dplyr::filter(.data$item_id %in% items_with_interactions$item_id) %>%
      dplyr::left_join(profs %>%
                         dplyr::rowwise() %>%
                         dplyr::mutate(total =
                                         sum(dplyr::c_across(where(is.numeric)))) %>%
                         tibble::rowid_to_column("profile") %>%
                         dplyr::filter(.data$total > 1) %>%
                         dplyr::select(-.data$total) %>%
                         tidyr::pivot_longer(cols = c(-.data$profile),
                                             names_to = "att",
                                             values_to = "mastered") %>%
                         dplyr::mutate(att =
                                         stringr::str_replace(.data$att,
                                                              "att_",
                                                              "mastered_")),
                       by = "profile") %>%
      dplyr::filter(!is.na(.data$att)) %>%
      dplyr::mutate(mastered_att =
                      as.numeric(stringr::str_remove(.data$att,
                                                     "mastered_"))) %>%
      dplyr::select(-.data$att) %>%
      dplyr::left_join(q_matrix %>%
                         tibble::rowid_to_column("item_id") %>%
                         tidyr::pivot_longer(cols = c(-.data$item_id),
                                             names_to = "att",
                                             values_to = "measured") %>%
                         dplyr::mutate(measured_att =
                                         as.numeric(stringr::str_remove(.data$att,
                                                                        "att_"))) %>%
                         dplyr::select(-.data$att),
                       by = "item_id") %>%
      dplyr::filter(.data$mastered_att == .data$measured_att) %>%
      dplyr::mutate(master = as.numeric(.data$mastered >= .data$measured)) %>%
      dplyr::group_by(.data$profile, .data$item_id) %>%
      dplyr::mutate(master = mean(.data$master)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$mastered, -.data$mastered_att) %>%
      dplyr::mutate(measured = .data$measured * .data$measured_att,
             measured_att =
               stringr::str_c("att_", as.character(.data$measured_att))) %>%
      dplyr::filter(.data$measured != 0) %>%
      dplyr::group_by(.data$profile, .data$item_id) %>%
      dplyr::mutate(meas =
                      stringr::str_c("att_",
                                     as.character(dplyr::row_number()))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$measured_att) %>%
      tidyr::pivot_wider(names_from = "meas", values_from = "measured") %>%
      dplyr::mutate(param = dplyr::case_when(.data$master < 1 ~ NA_character_,
                                             .data$master == 1 ~
                                               as.character(glue::glue("l{item_id}_2{att_1}{att_2}")))) %>%
      dplyr::select(.data$profile, .data$item_id, .data$param)
  } else {
    profile_item_interactions <-
      tibble::tibble(profile = rep(1:(2^ncol(q_matrix)),
                                   each = nrow(q_matrix)),
                     item_id = rep(seq_len(nrow(q_matrix)),
                                   times = (2^ncol(q_matrix)))) %>%
      dplyr::mutate(param = NA_character_)
  }

  pi_mat <- tibble::tibble(profile = rep(1:(2^ncol(q_matrix)),
                                         each = nrow(q_matrix)),
                           item_id = rep(seq_len(nrow(q_matrix)),
                                         times = (2^ncol(q_matrix)))) %>%
    dplyr::left_join(profs %>%
                       tibble::rowid_to_column("profile") %>%
                       tidyr::pivot_longer(cols = c(-.data$profile),
                                           names_to = "att_mastered",
                                           values_to = "mastered"),
                     by = "profile") %>%
    dplyr::left_join(q_matrix %>%
                       tibble::rowid_to_column("item_id") %>%
                       tidyr::pivot_longer(cols = c(-.data$item_id),
                                           names_to = "att_measured",
                                           values_to = "measured"),
                     by = "item_id") %>%
    dplyr::filter(.data$att_mastered == .data$att_measured) %>%
    dplyr::mutate(int0 = glue::glue("l_0"),
                  need_param = .data$mastered * .data$measured,
                  attribute =
                    as.numeric(stringr::str_remove(.data$att_measured, "att_")),
                  mef = dplyr::case_when(.data$need_param == 0 ~ NA_character_,
                                         .data$need_param > 0 ~
                                           as.character(glue::glue("l_1")))) %>%
    dplyr::select(-.data$att_measured, -.data$attribute, -.data$measured,
                  -.data$mastered, -.data$need_param) %>%
    tidyr::pivot_wider(names_from = "att_mastered", values_from = "mef") %>%
    dplyr::left_join(profile_item_interactions %>%
                       dplyr::rename(int2 = .data$param),
                     by = c("profile", "item_id")) %>%
    tidyr::unite(col = "param", c(-.data$profile, -.data$item_id), sep = "+",
                 na.rm = T) %>%
    dplyr::mutate(stan_pi =
             as.character(glue::glue("pi[{item_id},{profile}] = inv_logit({param});")))

  stan_data <-
    glue::glue("data {{",
               "  int<lower=1> I;                      // number of items",
               "  int<lower=1> J;                      // number of respondents",
               "  int<lower=1> N;                      // number of observations",
               "  int<lower=1> C;                      // number of classes",
               "  int<lower=1> A;                      // number of attributes",
               "  array[N, 2] int<lower=1,upper=I> ii; // item for obs n",
               "  array[N, 2] int<lower=0,upper=1> y;  // score for obs n",
               "  array[J, 2] int<lower=1,upper=N> s;  // starting row for j",
               "  array[J, 2] int<lower=1,upper=I> l;  // number of items for j",
               "  matrix[C,A] Alpha;                   // attribute pattern for each C",
               "}}", .sep = "\n")

  if (all(int2 == "")) {
    stan_parameters <-
      glue::glue("parameters {{",
                 "  array[C] simplex[C] tau;",
                 "  simplex[C] Vc;",
                 glue::glue_collapse(glue::glue("  {int0}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef}"), "\n"),
                 "}}", .sep = "\n")
  } else {
    stan_parameters <-
      glue::glue("parameters {{",
                 "  array[C] simplex[C] tau;",
                 "  simplex[C] Vc;",
                 glue::glue_collapse(glue::glue("  {int0}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef}"), "\n"),
                 glue::glue_collapse(glue::glue("  {int2}"), "\n"),
                 "}}", .sep = "\n")
  }

  stan_transformed_parameters <-
    glue::glue("transformed parameters {{",
               "  matrix[I,C] pi;",
               "",
               glue::glue_collapse(glue::glue("  {pi_mat$stan_pi}"), "\n"),
               "}}", .sep = "\n")

  if (all(int2_priors == "")) {
    stan_model <-
      glue::glue("model {{",
                 "  array[C, C] real ps;",
                 "",
                 "  // Priors",
                 glue::glue_collapse(glue::glue("  {int0_priors}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef_priors}"), "\n"),
                 "",
                 "  // Likelihood",
                 "  for (j in 1:J) {{",
                 "    vector[C] tmp;",
                 "    for (c1 in 1:C) {{",
                 "      for (c2 in 1:C) {{",
                 "        array[l[j, 1]] real log_items;",
                 "        for (m in 1:l[j, 1]) {{",
                 "          int i = ii[s[j, 1] + m - 1, 1];",
                 "          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);",
                 "        }}",
                 "        ps[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);",
                 "      }}",
                 "      tmp[c1] = log_sum_exp(ps[c1,]);",
                 "    }}",
                 "    target += log_sum_exp(tmp);",
                 "  }}",
                 "}}", .sep = "\n")
  } else {
    stan_model <-
      glue::glue("model {{",
                 "  array[C, C] real ps;",
                 "",
                 "  // Priors",
                 glue::glue_collapse(glue::glue("  {int0_priors}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef_priors}"), "\n"),
                 glue::glue_collapse(glue::glue("  {int2_priors}"), "\n"),
                 "",
                 "  // Likelihood",
                 "  for (j in 1:J) {{",
                 "    vector[C] tmp;",
                 "    for (c1 in 1:C) {{",
                 "      for (c2 in 1:C) {{",
                 "        array[l[j, 1]] real log_items;",
                 "        for (m in 1:l[j, 1]) {{",
                 "          int i = ii[s[j, 1] + m - 1, 1];",
                 "          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);",
                 "        }}",
                 "        ps[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);",
                 "      }}",
                 "      tmp[c1] = log_sum_exp(ps[c1,]);",
                 "    }}",
                 "    target += log_sum_exp(tmp);",
                 "  }}",
                 "}}", .sep = "\n")
  }

  stan_generated_quantities <-
    glue::glue("generated quantities {{",
               "  vector[J] log_lik;",
               "  array[J] matrix[C, C] prob_transition_class;",
               "  array[J] matrix[A, 2] prob_resp_attr;",
               "",
               "  // Likelihood",
               "  for (j in 1:J) {{",
               "    vector[C] tmp;",
               "    array[C, C] real ps;",
               "    for (c1 in 1:C) {{",
               "      for (c2 in 1:C) {{",
               "        array[l[j, 1]] real log_items;",
               "        for (m in 1:l[j, 1]) {{",
               "          int i = ii[s[j, 1] + m - 1, 1];",
               "          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);",
               "        }}",
               "        ps[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);",
               "      }}",
               "      tmp[c1] = log_sum_exp(ps[c1,]);",
               "    }}",
               "    log_lik[j] = log_sum_exp(tmp);",
               "  }}",
               "",
               "  // latent class probabilities",
               "  for (j in 1:J) {{",
               "    vector[C] tmp;",
               "    matrix[C, C] prob_joint;",
               "    for (c1 in 1:C) {{",
               "      for (c2 in 1:C) {{",
               "        array[l[j, 1]] real log_items;",
               "        for (m in 1:l[j, 1]) {{",
               "          int i = ii[s[j, 1] + m - 1, 1];",
               "          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);",
               "        }}",
               "        prob_joint[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);",
               "      }}",
               "    }}",
               "    prob_transition_class[j] = exp(prob_joint) / sum(exp(prob_joint));",
               "  }}",
               "",
               "  for (j in 1:J) {{",
               "    for (a in 1:A) {{",
               "      vector[C] prob_attr_class_t1;",
               "      vector[C] prob_attr_class_t2;",
               "      for (c in 1:C) {{",
               "        prob_attr_class_t1[c] = sum(prob_transition_class[j,c,]) * Alpha[c,a];",
               "        prob_attr_class_t2[c] = sum(prob_transition_class[j,,c]) * Alpha[c,a];",
               "      }}",
               "      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);",
               "      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);",
               "    }}",
               "  }}",
               "}}", .sep = "\n")

  stan_code <- list(data = stan_data,
                    parameters = stan_parameters,
                    transformed_parameters = stan_transformed_parameters,
                    model = stan_model,
                    generated_quantities = stan_generated_quantities)

  return(stan_code)
}
