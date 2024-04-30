#' Creating Multi-Threaded TDCM Stan Code
#'
#' Automating the creation of multi-threaded Stan code for a TDCM.
#'
#' @param q_matrix A tibble containing the assessment Q-matrix.
#'
#' @return `stan_code` A list containing the text for the Stan code blocks.
#'
#' @export
#'
#' @examples
#' qmatrix = tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 1),
#'                          att_2 = c(0, 1, 0, 1, 1, 1))
#' create_threaded_stan_tdcm(q_matrix = qmatrix)
create_threaded_stan_tdcm <- function(q_matrix) {
  profs <- bin_profile(ncol(q_matrix))

  colnames(q_matrix) <- glue::glue("att_{1:ncol(q_matrix)}")

  int0 <- glue::glue("real l{1:nrow(q_matrix)}_0;")
  int0_priors <- glue::glue("l{1:nrow(q_matrix)}_0 ~ normal(0, 2);")

  mef <- q_matrix %>%
    tibble::rowid_to_column("item_id") %>%
    tidyr::pivot_longer(cols = c(-"item_id"), names_to = "attr",
                        values_to = "meas") %>%
    dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                        "att_"))) %>%
    dplyr::filter(.data$meas == 1) %>%
    dplyr::select(-"meas") %>%
    dplyr::mutate(param = glue::glue("real<lower=0> l{item_id}_1{attr};")) %>%
    dplyr::pull(.data$param)
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
      tidyr::pivot_longer(cols = c(-"item_id"), names_to = "attr",
                          values_to = "meas") %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_num = dplyr::row_number(),
                    att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                               .data$att_num == 2 ~ "att2")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(attr =
                      as.numeric(stringr::str_remove(.data$attr, "att_"))) %>%
      dplyr::select(-"meas") %>%
      tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
      dplyr::mutate(param = glue::glue("real<lower=-1 * ",
                                       "fmin(l{item_id}_1{att1}, ",
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
      dplyr::mutate(attr =
                      as.numeric(stringr::str_remove(.data$attr, "att_"))) %>%
      dplyr::select(-"meas") %>%
      tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
      dplyr::mutate(param = glue::glue("l{item_id}_2{att1}{att2} ~ ",
                                       "normal(0, 2);")) %>%
      dplyr::pull(.data$param)
  }

  multi_item_q_matrix <- q_matrix %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(dplyr::c_across(where(is.numeric)))) %>%
    tibble::rowid_to_column("item_id") %>%
    dplyr::filter(.data$total == 2) %>%
    dplyr::select(-"total")

  if (nrow(multi_item_q_matrix) > 0) {
    items_with_interactions <- multi_item_q_matrix %>%
      tidyr::pivot_longer(cols = c(-"item_id"), names_to = "att",
                          values_to = "meas") %>%
      dplyr::mutate(meas_att =
                      as.numeric(stringr::str_remove(.data$att, "att_"))) %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_row = dplyr::row_number(),
                    att_row = stringr::str_c("att_",
                                             as.character(.data$att_row))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"att", -"meas") %>%
      tidyr::pivot_wider(names_from = "att_row", values_from = "meas_att") %>%
      dplyr::mutate(param =
                      as.character(glue::glue("l{item_id}_",
                                              "2{att_1}{att_2}"))) %>%
      dplyr::select("item_id", "param")

    profile_item_interactions <- tibble::tibble(profile =
                                                  rep(1:(2^ncol(q_matrix)),
                                                      each = nrow(q_matrix)),
                                                item_id =
                                                  rep(seq_len(nrow(q_matrix)),
                                                      times =
                                                      (2^ncol(q_matrix)))) %>%
      dplyr::filter(.data$item_id %in% items_with_interactions$item_id) %>%
      dplyr::left_join(profs %>%
                         dplyr::rowwise() %>%
                         dplyr::mutate(total =
                                         sum(dplyr::c_across(
                                           where(is.numeric)))) %>%
                         tibble::rowid_to_column("profile") %>%
                         dplyr::filter(.data$total > 1) %>%
                         dplyr::select(-"total") %>%
                         tidyr::pivot_longer(cols = c(-"profile"),
                                             names_to = "att",
                                             values_to = "mastered") %>%
                         dplyr::mutate(att =
                                         stringr::str_replace(.data$att,
                                                              "att_",
                                                              "mastered_")),
                       by = "profile", relationship = "many-to-many") %>%
      dplyr::filter(!is.na(.data$att)) %>%
      dplyr::mutate(mastered_att =
                      as.numeric(stringr::str_remove(.data$att,
                                                     "mastered_"))) %>%
      dplyr::select(-"att") %>%
      dplyr::left_join(q_matrix %>%
                         tibble::rowid_to_column("item_id") %>%
                         tidyr::pivot_longer(cols = c(-"item_id"),
                                             names_to = "att",
                                             values_to = "measured") %>%
                         dplyr::mutate(measured_att =
                                         as.numeric(stringr::str_remove(
                                           .data$att, "att_"))) %>%
                         dplyr::select(-"att"),
                       by = "item_id", relationship = "many-to-many") %>%
      dplyr::filter(.data$mastered_att == .data$measured_att) %>%
      dplyr::mutate(master = as.numeric(.data$mastered >= .data$measured)) %>%
      dplyr::group_by(.data$profile, .data$item_id) %>%
      dplyr::mutate(master = mean(.data$master)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"mastered", -"mastered_att") %>%
      dplyr::mutate(measured = .data$measured * .data$measured_att,
                    measured_att = stringr::str_c("att_",
                                                  as.character(
                                                    .data$measured_att))) %>%
      dplyr::filter(.data$measured != 0) %>%
      dplyr::group_by(.data$profile, .data$item_id) %>%
      dplyr::mutate(meas =
                      stringr::str_c("att_",
                                     as.character(dplyr::row_number()))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"measured_att") %>%
      tidyr::pivot_wider(names_from = "meas", values_from = "measured") %>%
      dplyr::mutate(param =
                      dplyr::case_when(
                        .data$master < 1 ~ NA_character_,
                        .data$master == 1 ~
                          as.character(glue::glue("l{item_id}_2",
                                                  "{att_1}{att_2}")))) %>%
      dplyr::select("profile", "item_id", "param")
  } else {
    profile_item_interactions <-
      tibble::tibble(profile = rep(1:(2^ncol(q_matrix)), each = nrow(q_matrix)),
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
                  attribute =
                    as.numeric(stringr::str_remove(.data$att_measured, "att_")),
                  mef = dplyr::case_when(.data$need_param == 0 ~ NA_character_,
                                         .data$need_param > 0 ~
                                           as.character(
                                             glue::glue(
                                               "l{item_id}_1{attribute}")))) %>%
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
               "    array[iis] int y1 = xi[1:iis];\n",
               "    array[iis] int ii1 = xi[(iis + 1):(2 * iis)];\n",
               "    array[iis] int jj1 = xi[((2 * iis) + 1):(3 * iis)];\n",
               "    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ",
               "ys)];\n",
               "    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + ",
               "(2 * ys))];\n",
               "\n",
               "    array[iis] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (2 * ys))];\n",
               "    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (3 * ys))];\n",
               "    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * ",
               "iis) + (4 * ys))];\n",
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
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
               "+ (C * C))];\n",
               "\n",
               "    matrix[C, C] ps;\n",
               "    matrix[C, C] tau_c;\n",
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
               "    for(c1 in 1:C) {{\n",
               "      for(c2 in 1:C) {{\n",
               "        int cc = c1 + ((c2 - 1) * C);\n",
               "        tau_c[c1, c2] = tauc[cc];\n",
               "      }}\n",
               "    }}\n",
               "\n",
               "    // Likelihood\n",
               "    for (j in 1:J) {{\n",
               "      vector[C] tmp;\n",
               "      for (c1 in 1:C) {{\n",
               "        for (c2 in 1:C) {{\n",
               "          array[l1[j]] real log_items;\n",
               "          for (m in 1:l1[j]) {{\n",
               "            int i = ii1[s1[j] + m - 1];\n",
               "            log_items[m] = y1[s1[j] + m - 1] * ",
               "log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - ",
               "pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - ",
               "y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);\n",
               "          }}\n",
               "          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + ",
               "sum(log_items);\n",
               "        }}\n",
               "        tmp[c1] = log_sum_exp(ps[c1,]);\n",
               "      }}\n",
               "      person += log_sum_exp(tmp);\n",
               "    }}\n",
               "\n",
               "    return [person]';\n",
               "  }}\n",
               "\n",
               "  vector person_loglik(vector beta, vector theta, array[] ",
               "real xr, array[] int xi) {{\n",
               "    int Z = num_elements(xi);\n",
               "    int ys = xi[Z - 1];\n",
               "    int iis = xi[Z];\n",
               "\n",
               "    array[iis] int y1 = xi[1:iis];\n",
               "    array[iis] int ii1 = xi[(iis + 1):(2 * iis)];\n",
               "    array[iis] int jj1 = xi[((2 * iis) + 1):(3 * iis)];\n",
               "    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ys)];\n",
               "    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + ",
               "(2 * ys))];\n",
               "\n",
               "    array[iis] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (2 * ys))];\n",
               "    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (3 * ys))];\n",
               "    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * ",
               "iis) + (4 * ys))];\n",
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
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
               "+ (C * C))];\n",
               "\n",
               "    matrix[C, C] ps;\n",
               "    matrix[C, C] tau_c;\n",
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
               "    for(c1 in 1:C) {{\n",
               "      for(c2 in 1:C) {{\n",
               "        int cc = c1 + ((c2 - 1) * C);\n",
               "        tau_c[c1, c2] = tauc[cc];\n",
               "      }}\n",
               "    }}\n",
               "\n",
               "    // Likelihood\n",
               "    for (j in 1:J) {{\n",
               "      vector[C] tmp;\n",
               "      for (c1 in 1:C) {{\n",
               "        for (c2 in 1:C) {{\n",
               "          array[l1[j]] real log_items;\n",
               "          for (m in 1:l1[j]) {{\n",
               "            int i = ii1[s1[j] + m - 1];\n",
               "            log_items[m] = y1[s1[j] + m - 1] * ",
               "log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - ",
               "pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - ",
               "y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);\n",
               "          }}\n",
               "          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + ",
               "sum(log_items);\n",
               "        }}\n",
               "        tmp[c1] = log_sum_exp(ps[c1,]);\n",
               "      }}\n",
               "      person[j] = log_sum_exp(tmp);\n",
               "    }}\n",
               "\n",
               "    return person;\n",
               "  }}\n",
               "\n",
               "  vector resp_transition(vector beta, vector theta, array[] ",
               "real xr, array[] int xi) {{\n",
               "    int Z = num_elements(xi);\n",
               "    int ys = xi[Z - 1];\n",
               "    int iis = xi[Z];\n",
               "\n",
               "    array[iis] int y1 = xi[1:iis];\n",
               "    array[iis] int ii1 = xi[(iis + 1):(2 * iis)];\n",
               "    array[iis] int jj1 = xi[((2 * iis) + 1):(3 * iis)];\n",
               "    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ys)];\n",
               "    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + ",
               "(2 * ys))];\n",
               "\n",
               "    array[iis] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * ",
               "iis) + (2 * ys))];\n",
               "    array[iis] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (2 * ys))];\n",
               "    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * ",
               "iis) + (3 * ys))];\n",
               "    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * ",
               "iis) + (4 * ys))];\n",
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
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) ",
               "+ (C * C))];\n",
               "\n",
               "    matrix[C, C] ps;\n",
               "    matrix[C, C] tau_c;\n",
               "    array[J] matrix[C, C] prob_transition_class;\n",
               "\n",
               "    vector[J * C * C] person;\n",
               "\n",
               "    matrix[I, C] pi_c;\n",
               "    for(c in 1:C) {{\n",
               "      for(i in 1:I) {{\n",
               "        int ic = i + ((c-1) * I);\n",
               "        pi_c[i, c] = pic[ic];\n",
               "      }}\n",
               "    }}\n",
               "\n",
               "    for(c1 in 1:C) {{\n",
               "      for(c2 in 1:C) {{\n",
               "        int cc = c1 + ((c2 - 1) * C);\n",
               "        tau_c[c1, c2] = tauc[cc];\n",
               "      }}\n",
               "    }}\n",
               "\n",
               "    // latent class probabilities\n",
               "    for (j in 1:J) {{\n",
               "      vector[C] tmp;\n",
               "      matrix[C, C] prob_joint;\n",
               "      for (c1 in 1:C) {{\n",
               "        for (c2 in 1:C) {{\n",
               "          array[l1[j]] real log_items;\n",
               "          for (m in 1:l1[j]) {{\n",
               "            int i = ii1[s1[j] + m - 1];\n",
               "            log_items[m] = y1[s1[j] + m - 1] * ",
               "log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - ",
               "pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - ",
               "y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);\n",
               "          }}\n",
               "          prob_joint[c1, c2] = log(Vc[c1]) + ",
               "log(tau_c[c1, c2]) + sum(log_items);\n",
               "        }}\n",
               "      }}\n",
               "      prob_transition_class[j] = exp(prob_joint) / ",
               "sum(exp(prob_joint));\n",
               "\n",
               "      for (c1 in 1:C) {{\n",
               "        for (c2 in 1:C) {{\n",
               "          person[((c1 - 1) * 8) + c2 + ((j - 1) * C * C)] ",
               "= prob_transition_class[j,c1,c2];\n",
               "        }}\n",
               "      }}\n",
               "    }}\n",
               "    return person;\n",
               "  }}\n",
               "}}\n", .sep = "")

  stan_data <-
    glue::glue("data {{",
               "  int<lower=1> I;",
               "  int<lower=1> J;",
               "  int<lower=1> N;",
               "  int<lower=1> C;",
               "  int<lower=1> A;",
               "  array[N, 2] int<lower=1,upper=I> ii;",
               "  array[N, 2] int<lower=1,upper=J> jj;",
               "  array[N, 2] int<lower=0,upper=1> y;",
               "  array[J, 2] int<lower=1,upper=N> s;",
               "  array[J, 2] int<lower=1,upper=I> l;",
               "  matrix[C,A] Alpha;",
               "  int<lower=1> n_shards;",
               "}}", .sep = "\n")

  stan_transformed_data <-
    glue::glue("transformed data {{\n",
               "  int ys = num_elements(s) / 2 / n_shards;\n",
               "  int iis = num_elements(ii) / 2 / n_shards;\n",
               "\n",
               "  int M = iis;\n",
               "\n",
               "  array[n_shards, (4 * ys) + (6 * iis) + 8] int xi;\n",
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
               "    xi[i, 1:iis] = y[iilower:iiupper, 1];\n",
               "    xi[i, (iis + 1):(iis + iis)] = ii[iilower:iiupper, 1];\n",
               "    xi[i, ((2 * iis) + 1):((2 * iis) + iis)] = ",
               "jj[iilower:iiupper, 1];\n",
               "    xi[i, ((3 * iis) + 1):((3 * iis) + ys)] = s[1:ys, 1];\n",
               "    xi[i, ((3 * iis) + ys + 1):((3 * iis) + (2 * ys))] = ",
               "l[ylower:yupper, 1];\n",
               "    xi[i, ((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ",
               "ys))] = y[iilower:iiupper, 2];\n",
               "    xi[i, ((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ",
               "ys))] = ii[iilower:iiupper, 2];\n",
               "    xi[i, ((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ",
               "ys))] = jj[iilower:iiupper, 2];\n",
               "    xi[i, ((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ",
               "ys))] = s[1:ys, 2];\n",
               "    xi[i, ((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ",
               "ys))] = l[ylower:yupper, 2];\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 1)] = I;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 2)] = N / n_shards;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 3)] = C;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 4)] = A;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 5)] = J / n_shards;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 6)] = iis;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 7)] = ys;\n",
               "    xi[i, ((6 * iis) + (4 * ys) + 8)] = iis;\n",
               "  }}\n",
               "}}\n", .sep = "")

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
    glue::glue("transformed parameters {{\n",
               "  matrix[I,C] pi;\n",
               "\n",
               glue::glue_collapse(glue::glue("  {pi_mat$stan_pi}"), "\n"),
               "\n",
               "\n",
               "  array[I * C] real pic;\n",
               "  for(c in 1:C) {{\n",
               "    for(i in 1:I) {{\n",
               "      int ic = i + ((c - 1) * I);\n",
               "      pic[ic] = pi[i, c];\n",
               "    }}\n",
               "  }}\n",
               "\n",
               "  array[C * C] real tauc;\n",
               "  for(c1 in 1:C) {{\n",
               "    for(c2 in 1:C) {{\n",
               "      int cc = c2 + ((c1 - 1) * C);\n",
               "      tauc[cc] = tau[c2, c1];\n",
               "    }}\n",
               "  }}\n",
               "\n",
               "  // a set of shared parameters\n",
               "  vector[C + (I * C) + (C * C)] beta;\n",
               "  beta[1:C] = Vc[1:C];\n",
               "  beta[(C + 1):(C + (I * C))] = to_vector(pic[1:(I * C)]);\n",
               "  beta[(C + (I * C) + 1):(C + (I * C) + (C * C))] = ",
               "to_vector(tauc[1:(C * C)]);\n",
               "}}\n", .sep = "")

  if (all(int2_priors == "")) {
    stan_model <-
      glue::glue("model {{",
                 "  array[C, C] real ps;",
                 "",
                 "  // Priors",
                 glue::glue_collapse(glue::glue("  {int0_priors}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef_priors}"), "\n"),
                 "",
                 "  target += sum(map_rect(sum_probs, beta, theta, xr, xi));",
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
                 "  target += sum(map_rect(sum_probs, beta, theta, xr, xi));",
                 "}}", .sep = "\n")
  }

  stan_generated_quantities <-
    glue::glue("generated quantities {{\n",
               "  vector[J] log_lik;\n",
               "  array[J] matrix[C, C] format_prob_transition_class;\n",
               "  vector[J*C*C] prob_transition_class;\n",
               "  array[J] matrix[A, 2] prob_resp_attr;\n",
               "\n",
               "  log_lik = map_rect(person_loglik, beta, theta, xr, xi);\n",
               "\n",
               "  prob_transition_class = map_rect(resp_transition, beta, ",
               "theta, xr, xi);\n",
               "\n",
               "  for (j in 1:J) {{\n",
               "    for (c1 in 1:C) {{\n",
               "      for (c2 in 1:C) {{\n",
               "        int iter = ((c1 - 1) * 8) + c2 + ((j - 1) * C * C);\n",
               "        format_prob_transition_class[j,c1,c2] = ",
               "prob_transition_class[iter];\n",
               "      }}\n",
               "    }}\n",
               "  }}\n",
               "\n",
               "  for (j in 1:J) {{\n",
               "    for (a in 1:A) {{\n",
               "      vector[C] prob_attr_class_t1;\n",
               "      vector[C] prob_attr_class_t2;\n",
               "      for (c in 1:C) {{\n",
               "        prob_attr_class_t1[c] = ",
               "sum(format_prob_transition_class[j,c,]) * Alpha[c,a];\n",
               "        prob_attr_class_t2[c] = ",
               "sum(format_prob_transition_class[j,,c]) * Alpha[c,a];\n",
               "      }}\n",
               "      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);\n",
               "      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);\n",
               "    }}\n",
               "  }}\n",
               "}}\n", .sep = "")

  stan_code <- list(functions = stan_functions,
                    data = stan_data,
                    transformed_data = stan_transformed_data,
                    parameters = stan_parameters,
                    transformed_parameters = stan_transformed_parameters,
                    model = stan_model,
                    generated_quantities = stan_generated_quantities)

  return(stan_code)
}
