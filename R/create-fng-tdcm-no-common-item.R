#' Creating Fungible TDCM with No Common Items Stan Code
#'
#' Automating the creation of fungible Stan code for a TDCM when there are no
#' common items.
#'
#' @param q_matrix A tibble containing the assessment Q-matrix.
#'
#' @return `stan_code` A list containing the text for the Stan code blocks.
#'
#' @export
#'
#' @examples
#' qmatrix = tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 0),
#'                          att_2 = c(0, 1, 0, 1, 0, 1))
#' create_fng_no_common_item_tdcm(q_matrix = qmatrix)
create_fng_no_common_item_tdcm <- function(q_matrix) {
  q_matrix <- check_multi_attribute(q_matrix)

  profs <- bin_profile(ncol(q_matrix))

  colnames(q_matrix) <- glue::glue("att_{1:ncol(q_matrix)}")

  int0 <- glue::glue("real l_0;")
  int0_priors <- glue::glue("l_0 ~ normal(0, 2);")

  mef <- glue::glue("real<lower=0> l_1;")
  mef_priors <- glue::glue("l_1 ~ lognormal(0, 1);")

  profile_item_interactions <- tibble::tibble(profile =
                                                rep(1:(2^ncol(q_matrix)),
                                                    each = nrow(q_matrix)),
                                              item_id =
                                                rep(seq_len(nrow(q_matrix)),
                                                    times =
                                                    (2^ncol(q_matrix)))) %>%
    dplyr::mutate(param = NA_character_)

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
    dplyr::mutate(int0 = glue::glue("l_0"), need_param = .data$mastered *
                    .data$measured,
                  attribute = as.numeric(stringr::str_remove(.data$att_measured,
                                                             "att_")),
                  mef = dplyr::case_when(.data$need_param ==
                                           0 ~ NA_character_,
                                         .data$need_param > 0 ~
                                           as.character(glue::glue("l_1")))) %>%
    dplyr::select(-"att_measured", -"attribute",
                  -"measured", -"mastered", -"need_param") %>%
    tidyr::pivot_wider(names_from = "att_mastered", values_from = "mef") %>%
    dplyr::left_join(profile_item_interactions %>%
                       dplyr::rename(int2 = "param"),
                     by = c("profile", "item_id"),
                     relationship = "many-to-many") %>%
    tidyr::unite(col = "param", c(-"profile", -"item_id"), sep = "+",
                 na.rm = TRUE) %>%
    dplyr::mutate(stan_pi = as.character(glue::glue("pi[{item_id},{profile}] =",
                                                    " inv_logit({param});")))
  stan_data <- glue::glue("data {{",
                          "  int<lower=1> I;",
                          "  int<lower=1> J;",
                          "  int<lower=1> N;",
                          "  int<lower=1> C;",
                          "  int<lower=1> A;",
                          "  array[N, 2] int<lower=1,upper=I> ii;",
                          "  array[N, 2] int<lower=0> y;",
                          "  array[J, 2] int<lower=1,upper=N> s;",
                          "  array[J, 2] int<lower=1,upper=I> l;",
                          "  matrix[C,A] Alpha;",
                          "}}", .sep = "\n")

  stan_parameters <- glue::glue("parameters {{",
                                "  array[C] simplex[C] tau;",
                                "  simplex[C] Vc;",
                                glue::glue_collapse(glue::glue("  {int0}"),
                                                    "\n"),
                                glue::glue_collapse(glue::glue("  {mef}"),
                                                    "\n"),
                                "}}", .sep = "\n")

  stan_transformed_parameters <-
    glue::glue("transformed parameters {{",
               "  matrix[I,C] pi;",
               "",
               glue::glue_collapse(glue::glue("  {pi_mat$stan_pi}"), "\n"),
               "}}", .sep = "\n")

  stan_model <- glue::glue("model {{\n",
                           "  array[C, C] real ps;\n",
                           "\n",
                           "  // Priors\n",
                           glue::glue_collapse(glue::glue("  {int0_priors}"),
                                               "\n"),
                           "\n",
                           glue::glue_collapse(glue::glue("  {mef_priors}"),
                                               "\n"),
                           "\n",
                           "\n",
                           "  // Likelihood\n",
                           "  for (j in 1:J) {{\n",
                           "    vector[C] tmp;\n",
                           "    for (c1 in 1:C) {{\n",
                           "      for (c2 in 1:C) {{\n",
                           "        array[l[j, 1]] real log_items;\n",
                           "        for (m in 1:l[j, 1]) {{\n",
                           "          int i = ii[s[j, 1] + m - 1, 1];\n",
                           "          real tmp1 = 0;\n",
                           "          real tmp2 = 0;\n",
                           "          if(y[s[j, 1] + m - 1, 1] != 9) ",
                           "{{tmp1 = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) +",
                           " (1 - y[s[j, 1] + m - 1, 1]) * ",
                           "log(1 - pi[i,c1]);}}\n",
                           "          if(y[s[j, 1] + m - 1, 2] != 9) ",
                           "{{tmp2 = y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) ",
                           "+ (1 - y[s[j, 1] + m - 1, 2]) * ",
                           "log(1 - pi[i,c2]);}}\n",
                           "          log_items[m] = tmp1 + tmp2;\n",
                           "        }}\n",
                           "        ps[c1, c2] = log(Vc[c1]) + ",
                           "log(tau[c1, c2]) + sum(log_items);\n",
                           "      }}\n",
                           "      tmp[c1] = log_sum_exp(ps[c1,]);\n",
                           "    }}\n",
                           "    target += log_sum_exp(tmp);\n",
                           "  }}\n",
                           "}}\n", .sep = "")

  stan_generated_quantities <- glue::glue("generated quantities {{\n",
                                          "  vector[J] log_lik;\n",
                                          "  array[J] matrix[C, C] ",
                                          "prob_transition_class;\n",
                                          "  array[J] matrix[A, 2] ",
                                          "prob_resp_attr;\n",
                                          "  array[N, 2] int<lower=0> y_rep;\n",
                                          "  array[J, 2] int j_class;\n",
                                          "  vector[C] tau_tmp;\n",
                                          "\n",
                                          "  // Likelihood\n",
                                          "  for (j in 1:J) {{\n",
                                          "    vector[C] tmp;\n",
                                          "    array[C, C] real ps;\n",
                                          "    for (c1 in 1:C) {{\n",
                                          "      for (c2 in 1:C) {{\n",
                                          "        array[l[j, 1]] real ",
                                          "log_items;\n",
                                          "        for (m in 1:l[j, 1]) {{\n",
                                          "          int i = ii[s[j, 1] + ",
                                          "m - 1, 1];\n",
                                          "          real tmp1 = 0;\n",
                                          "          real tmp2 = 0;\n",
                                          "          if(y[s[j, 1] + m - 1, 1] ",
                                          "!= 9) {{tmp1 = y[s[j, 1] + m - 1, ",
                                          "1] * log(pi[i,c1]) + (1 - ",
                                          "y[s[j, 1] + m - 1, 1]) * ",
                                          "log(1 - pi[i,c1]);}}\n",
                                          "          if(y[s[j, 1] + m - 1, 2] ",
                                          "!= 9) {{tmp2 = y[s[j, 1] + m - 1, ",
                                          "2] * log(pi[i,c2]) + (1 - ",
                                          "y[s[j, 1] + m - 1, 2]) * ",
                                          "log(1 - pi[i,c2]);}}\n",
                                          "          log_items[m] = tmp1 + ",
                                          "tmp2;\n",
                                          "        }}\n",
                                          "        ps[c1, c2] = log(Vc[c1]) + ",
                                          "log(tau[c1, c2]) + ",
                                          "sum(log_items);\n",
                                          "      }}\n",
                                          "      tmp[c1] = ",
                                          "log_sum_exp(ps[c1,]);\n",
                                          "    }}\n",
                                          "    log_lik[j] = ",
                                          "log_sum_exp(tmp);\n",
                                          "  }}\n",
                                          "\n",
                                          "  // latent class probabilities\n",
                                          "  for (j in 1:J) {{\n",
                                          "    vector[C] tmp;\n",
                                          "    matrix[C, C] prob_joint;\n",
                                          "    for (c1 in 1:C) {{\n",
                                          "      for (c2 in 1:C) {{\n",
                                          "        array[l[j, 1]] real ",
                                          "log_items;\n",
                                          "        for (m in 1:l[j, 1]) {{\n",
                                          "          int i = ii[s[j, 1] + m - ",
                                          "1, 1];\n",
                                          "          real tmp1 = 0;\n",
                                          "          real tmp2 = 0;\n",
                                          "          if(y[s[j, 1] + m - 1, 1] ",
                                          "!= 9) {{tmp1 = y[s[j, 1] + m - 1, ",
                                          "1] * log(pi[i,c1]) + (1 - y[s[j, ",
                                          "1] + m - 1, 1]) * log(1 - ",
                                          "pi[i,c1]);}}\n",
                                          "          if(y[s[j, 1] + m - 1, ",
                                          "2] != 9) {{tmp2 = y[s[j, 1] + m - ",
                                          "1, 2] * log(pi[i,c2]) + (1 - ",
                                          "y[s[j, 1] + m - 1, 2]) * log(1 - ",
                                          "pi[i,c2]);}}\n",
                                          "          log_items[m] = tmp1 + ",
                                          "tmp2;\n",
                                          "        }}\n",
                                          "        prob_joint[c1, c2] = ",
                                          "log(Vc[c1]) + log(tau[c1, c2]) + ",
                                          "sum(log_items);\n",
                                          "      }}\n",
                                          "    }}\n",
                                          "    prob_transition_class[j] = ",
                                          "exp(prob_joint) / ",
                                          "sum(exp(prob_joint));\n",
                                          "  }}\n",
                                          "\n",
                                          "  for (j in 1:J) {{\n",
                                          "    for (a in 1:A) {{\n",
                                          "      vector[C] ",
                                          "prob_attr_class_t1;\n",
                                          "      vector[C] ",
                                          "prob_attr_class_t2;\n",
                                          "      for (c in 1:C) {{\n",
                                          "        prob_attr_class_t1[c] = ",
                                          "sum(prob_transition_class[j,c,]) * ",
                                          "Alpha[c,a];\n",
                                          "        prob_attr_class_t2[c] = ",
                                          "sum(prob_transition_class[j,,c]) * ",
                                          "Alpha[c,a];\n",
                                          "      }}\n",
                                          "      prob_resp_attr[j,a,1] = ",
                                          "sum(prob_attr_class_t1);\n",
                                          "      prob_resp_attr[j,a,2] = ",
                                          "sum(prob_attr_class_t2);\n",
                                          "    }}\n",
                                          "  }}\n",
                                          "\n",
                                          "  for (j in 1:J) {{\n",
                                          "    vector[C] j_probs = Vc / sum(Vc);\n",
                                          "    j_class[j, 1] = categorical_rng(j_probs);\n",
                                          "    tau_tmp = tau[j_class[j, 1], ];\n",
                                          "    j_class[j, 2] = categorical_rng(j_probs);\n",
                                          "    for (t in 1:2) {{\n",
                                          "      for (m in 1:l[j, t]) {{\n",
                                          "        int i = ii[s[j, t] + m - 1, t];\n",
                                          "        y_rep[x[j, t] + m - 1, t] = bernoulli_rng(pi[i, j_class[j, t]]);\n",
                                          "      }}\n",
                                          "    }}\n",
                                          "  }}\n",
                                          "}}\n", .sep = "")

  stan_code <- list(data = stan_data, parameters = stan_parameters,
                    transformed_parameters = stan_transformed_parameters,
                    model = stan_model,
                    generated_quantities = stan_generated_quantities)

  return(stan_code)
}
