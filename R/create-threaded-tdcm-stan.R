#' Creating Multi-Threaded TDCM Stan Code
#'
#' Automating the creation of multi-threaded Stan code for a TDCM.
#'
#' @param q_matrix A tibble containing the assessment Q-matrix.
#' @param profs A tibble of the possible attribute mastery profiles.
#'
#' @return `stan_code` A list containing the text for the Stan code blocks.
#'
#' @export
#'
#' @examples
#' qmatrix = tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 1), att_2 = c(0, 1, 0, 1, 1, 1))
#' possible_profiles = tibble::tibble(att_1 = c(0, 1, 0, 1), att_2 = c(0, 1, 0, 1))
#' create_threaded_stan_tdcm(q_matrix = qmatrix, profs = possible_profiles)
create_threaded_stan_tdcm <- function(q_matrix, profs) {
  colnames(q_matrix) <- glue::glue("att_{1:ncol(q_matrix)}")

  int0 <- glue::glue("real l{1:nrow(q_matrix)}_0;")
  int0_priors <- glue::glue("l{1:nrow(q_matrix)}_0 ~ normal(0, 2);")

  mef <- q_matrix %>%
    tibble::rowid_to_column("item_id") %>%
    tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "attr",
                        values_to = "meas") %>%
    dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                        "att_"))) %>%
    dplyr::filter(.data$meas == 1) %>%
    dplyr::select(-.data$meas) %>%
    dplyr::mutate(param = glue::glue("real<lower=0> l{item_id}_1{attr};")) %>%
    dplyr::pull(.data$param)
  mef_priors <- q_matrix %>%
    tibble::rowid_to_column("item_id") %>%
    tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "attr",
                        values_to = "meas") %>%
    dplyr::mutate(attr = as.numeric(stringr::str_remove(.data$attr,
                                                        "att_"))) %>%
    dplyr::filter(.data$meas == 1) %>%
    dplyr::select(-.data$meas) %>%
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
      tidyr::pivot_longer(cols = c(-.data$item_id), names_to = "attr",
                          values_to = "meas") %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_num = dplyr::row_number(),
                    att_num = dplyr::case_when(.data$att_num == 1 ~ "att1",
                                               .data$att_num == 2 ~ "att2")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(attr =
                      as.numeric(stringr::str_remove(.data$attr, "att_"))) %>%
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
      dplyr::mutate(attr =
                      as.numeric(stringr::str_remove(.data$attr, "att_"))) %>%
      dplyr::select(-.data$meas) %>%
      tidyr::pivot_wider(names_from = "att_num", values_from = "attr") %>%
      dplyr::mutate(param = glue::glue("l{item_id}_2{att1}{att2} ~ normal(0, 2);")) %>%
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
      dplyr::mutate(meas_att =
                      as.numeric(stringr::str_remove(.data$att, "att_"))) %>%
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

    profile_item_interactions <- tibble::tibble(profile =
                                                  rep(1:(2^ncol(q_matrix)),
                                                      each = nrow(q_matrix)),
                                                item_id = rep(seq_len(nrow(q_matrix)),
                                                              times = (2^ncol(q_matrix)))) %>%
      dplyr::filter(.data$item_id %in% items_with_interactions$item_id) %>%
      dplyr::left_join(profs %>%
                         dplyr::rowwise() %>%
                         dplyr::mutate(total = sum(dplyr::c_across(where(is.numeric)))) %>%
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
                    measured_att = stringr::str_c("att_",
                                                  as.character(.data$measured_att))) %>%
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
    dplyr::mutate(int0 = glue::glue("l{item_id}_0"),
                  need_param = .data$mastered * .data$measured,
                  attribute =
                    as.numeric(stringr::str_remove(.data$att_measured, "att_")),
                  mef = dplyr::case_when(.data$need_param == 0 ~ NA_character_,
                                         .data$need_param > 0 ~
                                           as.character(glue::glue("l{item_id}_1{attribute}")))) %>%
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

  stan_functions <-
    glue::glue("functions {{",
               "  real minmax (real x) {{",
               "    if (x < .01) {{",
               "      return 0.01;",
               "    }}",
               "",
               "    if (x > 0.99) {{",
               "      return 0.99;",
               "    }}",
               "",
               "    return x;",
               "  }}",
               "",
               "  vector sum_probs(vector beta, vector theta, real[] xr, int[] xi) {{",
               "    int Z = num_elements(xi);",
               "    int ys = xi[Z - 1];",
               "    int iis = xi[Z];",
               "",
               "    int y1[iis] = xi[1:iis];",
               "    int ii1[iis] = xi[(iis + 1):(2 * iis)];",
               "    int jj1[iis] = xi[((2 * iis) + 1):(3 * iis)];",
               "    int s1[ys] = xi[((3 * iis) + 1):((3 * iis) + ys)];",
               "    int l1[ys] = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];",
               "",
               "    int y2[iis] = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];",
               "    int ii2[iis] = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];",
               "    int jj2[iis] = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];",
               "    int s2[ys] = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];",
               "    int l2[ys] = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];",
               "",
               "    int I = xi[Z - 7];",
               "    int N = xi[Z - 6];",
               "    int C = xi[Z - 5];",
               "    int A = xi[Z - 4];",
               "    int J = xi[Z - 3];",
               "    int M = xi[Z - 2];",
               "",
               "    vector[C] Vc = beta[1:C];",
               "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];",
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];",
               "",
               "    matrix[C, C] ps;",
               "    matrix[C, C] tau_c;",
               "    real person = 0;",
               "",
               "    matrix[I, C] pi_c;",
               "    for(c in 1:C) {{",
               "      for(i in 1:I) {{",
               "        int ic = i + ((c-1) * I);",
               "        pi_c[i, c] = pic[ic];",
               "      }}",
               "    }}",
               "",
               "    for(c1 in 1:C) {{",
               "      for(c2 in 1:C) {{",
               "        int cc = c1 + ((c2 - 1) * C);",
               "        tau_c[c1, c2] = tauc[cc];",
               "      }}",
               "    }}",
               "",
               "    // Likelihood",
               "    for (j in 1:J) {{",
               "      vector[C] tmp;",
               "      for (c1 in 1:C) {{",
               "        for (c2 in 1:C) {{",
               "          real log_items[l1[j]];",
               "          for (m in 1:l1[j]) {{",
               "            int i = ii1[s1[j] + m - 1];",
               "            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);",
               "          }}",
               "          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);",
               "        }}",
               "        tmp[c1] = log_sum_exp(ps[c1,]);",
               "      }}",
               "      person += log_sum_exp(tmp);",
               "    }}",
               "",
               "    return [person]';",
               "  }}",
               "",
               "  vector person_loglik(vector beta, vector theta, real[] xr, int[] xi) {{",
               "    int Z = num_elements(xi);",
               "    int ys = xi[Z - 1];",
               "    int iis = xi[Z];",
               "",
               "    int y1[iis] = xi[1:iis];",
               "    int ii1[iis] = xi[(iis + 1):(2 * iis)];",
               "    int jj1[iis] = xi[((2 * iis) + 1):(3 * iis)];",
               "    int s1[ys] = xi[((3 * iis) + 1):((3 * iis) + ys)];",
               "    int l1[ys] = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];",
               "",
               "    int y2[iis] = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];",
               "    int ii2[iis] = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];",
               "    int jj2[iis] = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];",
               "    int s2[ys] = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];",
               "    int l2[ys] = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];",
               "",
               "    int I = xi[Z - 7];",
               "    int N = xi[Z - 6];",
               "    int C = xi[Z - 5];",
               "    int A = xi[Z - 4];",
               "    int J = xi[Z - 3];",
               "    int M = xi[Z - 2];",
               "",
               "    vector[C] Vc = beta[1:C];",
               "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];",
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];",
               "",
               "    matrix[C, C] ps;",
               "    matrix[C, C] tau_c;",
               "    vector[J] person;",
               "",
               "    matrix[I, C] pi_c;",
               "    for(c in 1:C) {{",
               "      for(i in 1:I) {{",
               "        int ic = i + ((c-1) * I);",
               "        pi_c[i, c] = pic[ic];",
               "      }}",
               "    }}",
               "",
               "    for(c1 in 1:C) {{",
               "      for(c2 in 1:C) {{",
               "        int cc = c1 + ((c2 - 1) * C);",
               "        tau_c[c1, c2] = tauc[cc];",
               "      }}",
               "    }}",
               "",
               "    // Likelihood",
               "    for (j in 1:J) {{",
               "      vector[C] tmp;",
               "      for (c1 in 1:C) {{",
               "        for (c2 in 1:C) {{",
               "          real log_items[l1[j]];",
               "          for (m in 1:l1[j]) {{",
               "            int i = ii1[s1[j] + m - 1];",
               "            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);",
               "          }}",
               "          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);",
               "        }}",
               "        tmp[c1] = log_sum_exp(ps[c1,]);",
               "      }}",
               "      person[j] = log_sum_exp(tmp);",
               "    }}",
               "",
               "    return person;",
               "  }}",
               "",
               "  vector resp_transition(vector beta, vector theta, real[] xr, int[] xi) {{",
               "    int Z = num_elements(xi);",
               "    int ys = xi[Z - 1];",
               "    int iis = xi[Z];",
               "",
               "    int y1[iis] = xi[1:iis];",
               "    int ii1[iis] = xi[(iis + 1):(2 * iis)];",
               "    int jj1[iis] = xi[((2 * iis) + 1):(3 * iis)];",
               "    int s1[ys] = xi[((3 * iis) + 1):((3 * iis) + ys)];",
               "    int l1[ys] = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];",
               "",
               "    int y2[iis] = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];",
               "    int ii2[iis] = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];",
               "    int jj2[iis] = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];",
               "    int s2[ys] = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];",
               "    int l2[ys] = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];",
               "",
               "    int I = xi[Z - 7];",
               "    int N = xi[Z - 6];",
               "    int C = xi[Z - 5];",
               "    int A = xi[Z - 4];",
               "    int J = xi[Z - 3];",
               "    int M = xi[Z - 2];",
               "",
               "    vector[C] Vc = beta[1:C];",
               "    vector[I * C] pic = beta[(C + 1):(C + (I * C))];",
               "    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];",
               "",
               "    matrix[C, C] ps;",
               "    matrix[C, C] tau_c;",
               "    matrix[C, C] prob_transition_class[J];",
               "",
               "    vector[J * C * C] person;",
               "",
               "    matrix[I, C] pi_c;",
               "    for(c in 1:C) {{",
               "      for(i in 1:I) {{",
               "        int ic = i + ((c-1) * I);",
               "        pi_c[i, c] = pic[ic];",
               "      }}",
               "    }}",
               "",
               "    for(c1 in 1:C) {{",
               "      for(c2 in 1:C) {{",
               "        int cc = c1 + ((c2 - 1) * C);",
               "        tau_c[c1, c2] = tauc[cc];",
               "      }}",
               "    }}",
               "",
               "    // latent class probabilities",
               "    for (j in 1:J) {{",
               "      vector[C] tmp;",
               "      matrix[C, C] prob_joint;",
               "      for (c1 in 1:C) {{",
               "        for (c2 in 1:C) {{",
               "          real log_items[l1[j]];",
               "          for (m in 1:l1[j]) {{",
               "            int i = ii1[s1[j] + m - 1];",
               "            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);",
               "          }}",
               "          prob_joint[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);",
               "        }}",
               "      }}",
               "      prob_transition_class[j] = exp(prob_joint) / sum(exp(prob_joint));",
               "",
               "      for (c1 in 1:C) {{",
               "        for (c2 in 1:C) {{",
               "          person[((c1 - 1) * 8) + c2 + ((j - 1) * C * C)] = prob_transition_class[j,c1,c2];",
               "        }}",
               "      }}",
               "    }}",
               "    return person;",
               "  }}",
               "}}", .sep = "\n")

  stan_data <-
    glue::glue("data {{",
               "  int<lower=1> I;                    // number of items",
               "  int<lower=1> J;                    // number of respondents",
               "  int<lower=1> N;                    // number of observations",
               "  int<lower=1> C;                    // number of classes",
               "  int<lower=1> A;                    // number of attributes",
               "  int<lower=1,upper=I> ii[N, 2];     // item for obs n",
               "  int<lower=1,upper=J> jj[N, 2];     // respondent for obs n",
               "  int<lower=0,upper=1> y[N, 2];      // score for obs n",
               "  int<lower=1,upper=N> s[J, 2];      // starting row for j",
               "  int<lower=1,upper=I> l[J, 2];      // number of items for j",
               "  matrix[C,A] Alpha;                 // attribute pattern for each C",
               "  int<lower=1> n_shards;             // the number of shards to split the data into",
               "}}", .sep = "\n")

  stan_transformed_data <-
    glue::glue("transformed data {{",
               "  int ys = num_elements(s) / 2 / n_shards;",
               "  int iis = num_elements(ii) / 2 / n_shards;",
               "",
               "  int M = iis;",
               "",
               "  int xi[n_shards, (4 * ys) + (6 * iis) + 8];",
               "",
               "  // an empty set of per-shard parameters",
               "  vector[0] theta[n_shards];",
               "",
               "  real xr[n_shards,1];",
               "  for(kk in 1:n_shards) {{",
               "    xr[kk, 1] = 1.0;",
               "  }}",
               "",
               "  // split into shards",
               "  for (i in 1:n_shards) {{",
               "    int ylower;",
               "    int yupper;",
               "    int iilower;",
               "    int iiupper;",
               "",
               "    ylower = ((i - 1) * ys) + 1;",
               "    yupper = i * ys;",
               "    iilower = ((i - 1) * iis) + 1;",
               "    iiupper = i * iis;",
               "",
               "    xi[i, 1:iis] = y[iilower:iiupper, 1];",
               "    xi[i, (iis + 1):(iis + iis)] = ii[iilower:iiupper, 1];",
               "    xi[i, ((2 * iis) + 1):((2 * iis) + iis)] = jj[iilower:iiupper, 1];",
               "    xi[i, ((3 * iis) + 1):((3 * iis) + ys)] = s[1:ys, 1];",
               "    xi[i, ((3 * iis) + ys + 1):((3 * iis) + (2 * ys))] = l[ylower:yupper, 1];",
               "    xi[i, ((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))] = y[iilower:iiupper, 2];",
               "    xi[i, ((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))] = ii[iilower:iiupper, 2];",
               "    xi[i, ((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))] = jj[iilower:iiupper, 2];",
               "    xi[i, ((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))] = s[1:ys, 2];",
               "    xi[i, ((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))] = l[ylower:yupper, 2];",
               "    xi[i, ((6 * iis) + (4 * ys) + 1)] = I;",
               "    xi[i, ((6 * iis) + (4 * ys) + 2)] = N / n_shards;",
               "    xi[i, ((6 * iis) + (4 * ys) + 3)] = C;",
               "    xi[i, ((6 * iis) + (4 * ys) + 4)] = A;",
               "    xi[i, ((6 * iis) + (4 * ys) + 5)] = J / n_shards;",
               "    xi[i, ((6 * iis) + (4 * ys) + 6)] = iis;",
               "    xi[i, ((6 * iis) + (4 * ys) + 7)] = ys;",
               "    xi[i, ((6 * iis) + (4 * ys) + 8)] = iis;",
               "  }}",
               "}}", .sep = "\n")

  if (all(int2 == "")) {
    stan_parameters <-
      glue::glue("parameters {{",
                 "  simplex[C] tau[C];",
                 "  simplex[C] Vc;",
                 glue::glue_collapse(glue::glue("  {int0}"), "\n"),
                 glue::glue_collapse(glue::glue("  {mef}"), "\n"),
                 "}}", .sep = "\n")
  } else {
    stan_parameters <-
      glue::glue("parameters {{",
                 "  simplex[C] tau[C];",
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
               "",
               "  real pic[I * C];",
               "  for(c in 1:C) {{",
               "    for(i in 1:I) {{",
               "      int ic = i + ((c - 1) * I);",
               "      pic[ic] = pi[i, c];",
               "    }}",
               "  }}",
               "",
               "  real tauc[C * C];",
               "  for(c1 in 1:C) {{",
               "    for(c2 in 1:C) {{",
               "      int cc = c2 + ((c1 - 1) * C);",
               "      tauc[cc] = tau[c2, c1];",
               "    }}",
               "  }}",
               "",
               "  // a set of shared parameters",
               "  vector[C + (I * C) + (C * C)] beta;",
               "  beta[1:C] = Vc[1:C];",
               "  beta[(C + 1):(C + (I * C))] = to_vector(pic[1:(I * C)]);",
               "  beta[(C + (I * C) + 1):(C + (I * C) + (C * C))] = to_vector(tauc[1:(C * C)]);",
               "}}", .sep = "\n")

  if (all(int2_priors == "")) {
    stan_model <-
      glue::glue("model {{",
                 "  real ps[C, C];",
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
                 "  real ps[C, C];",
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
    glue::glue("generated quantities {{",
               "  vector[J] log_lik;",
               "  matrix[C, C] format_prob_transition_class[J];",
               "  vector[J*C*C] prob_transition_class;",
               "  matrix[A, 2] prob_resp_attr[J];",
               "",
               "  log_lik = map_rect(person_loglik, beta, theta, xr, xi);",
               "",
               "  prob_transition_class = map_rect(resp_transition, beta, theta, xr, xi);",
               "",
               "  for (j in 1:J) {{",
               "    for (c1 in 1:C) {{",
               "      for (c2 in 1:C) {{",
               "        int iter = ((c1 - 1) * 8) + c2 + ((j - 1) * C * C);",
               "        format_prob_transition_class[j,c1,c2] = prob_transition_class[iter];",
               "      }}",
               "    }}",
               "  }}",
               "",
               "  for (j in 1:J) {{",
               "    for (a in 1:A) {{",
               "      vector[C] prob_attr_class_t1;",
               "      vector[C] prob_attr_class_t2;",
               "      for (c in 1:C) {{",
               "        prob_attr_class_t1[c] = sum(format_prob_transition_class[j,c,]) * Alpha[c,a];",
               "        prob_attr_class_t2[c] = sum(format_prob_transition_class[j,,c]) * Alpha[c,a];",
               "      }}",
               "      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);",
               "      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);",
               "    }}",
               "  }}",
               "}}", .sep = "\n")

  stan_code <- list(functions = stan_functions,
                    data = stan_data,
                    transformed_data = stan_transformed_data,
                    parameters = stan_parameters,
                    transformed_parameters = stan_transformed_parameters,
                    model = stan_model,
                    generated_quantities = stan_generated_quantities)

  return(stan_code)
}
