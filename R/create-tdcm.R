#' Creating TDCM Stan Code
#'
#' Automating the creation of Stan code for a TDCM.
#'
#' @param q_matrix A tibble containing the assessment Q-matrix.
#' @param t An intenger containing the number of assessment points.
#' @param multithread A logical variable indicating whether to using
#' multithreading.
#' @param fungible A logical variable indicating whether to estimate a fungible
#' model.
#' @param repeated_items A logical variable indicating whether the same
#' items were administered at each assessment point.
#'
#' @return `stan_code` A list containing the text for the Stan code blocks.
#'
#' @export
#'
#' @examples
#' qmatrix = tibble::tibble(att_1 = c(1, 0, 1, 0, 1, 1),
#'                          att_2 = c(0, 1, 0, 1, 1, 1))
#' create_stan_tdcm(q_matrix = qmatrix, t = 2, multithread = FALSE,
#' fungible = FALSE)
create_stan_tdcm <- function(q_matrix, t = 2, multithread = TRUE,
                             fungible = FALSE, repeated_items = TRUE) {
  profs <- bin_profile(ncol(q_matrix))

  colnames(q_matrix) <- glue::glue("att_{1:ncol(q_matrix)}")

  priors <- calc_priors(q_matrix, fungible)
  int0_priors <- priors$intercept
  mef_priors <- priors$main_effect
  int2_priors <- priors$interactions

  multi_item_q_matrix <- q_matrix %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total =
                    sum(dplyr::c_across(where(is.numeric)))) %>%
    dplyr::ungroup() %>%
    tibble::rowid_to_column("item_id") %>%
    dplyr::filter(.data$total == 2) %>%
    dplyr::select(-"total")

  if (fungible) {
    int0 <- glue::glue("real l_0;")
    mef <- glue::glue("real<lower=0> l_1;")
    int2 <- ""
  } else {
    int0 <- glue::glue("real l{1:nrow(q_matrix)}_0;")

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

    if (nrow(multi_item_q_matrix) == 0) {
      int2 <- ""
    } else {
      int2 <- multi_item_q_matrix %>%
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
    }
  }

  if (nrow(multi_item_q_matrix) > 0) {
    items_with_interactions <- multi_item_q_matrix %>%
      tidyr::pivot_longer(cols = c(-"item_id"), names_to = "att",
                          values_to = "meas") %>%
      dplyr::mutate(meas_att = as.numeric(stringr::str_remove(.data$att,
                                                              "att_"))) %>%
      dplyr::filter(.data$meas == 1) %>%
      dplyr::group_by(.data$item_id) %>%
      dplyr::mutate(att_row = dplyr::row_number(),
                    att_row = stringr::str_c("att_",
                                             as.character(.data$att_row))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"att", -"meas") %>%
      tidyr::pivot_wider(names_from = "att_row", values_from = "meas_att") %>%
      dplyr::mutate(param =
                      as.character(glue::glue("l{item_id}_2",
                                              "{att_1}{att_2}"))) %>%
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
                         dplyr::mutate(total = sum(dplyr::c_across(
                           where(is.numeric)))) %>%
                         tibble::rowid_to_column("profile") %>%
                         dplyr::filter(.data$total > 1) %>%
                         dplyr::select(-"total") %>%
                         tidyr::pivot_longer(cols = c(-"profile"),
                                             names_to = "att",
                                             values_to = "mastered") %>%
                         dplyr::mutate(att =
                                         stringr::str_replace(.data$att, "att_",
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
                                         as.numeric(
                                           stringr::str_remove(
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
                    measured_att =
                      stringr::str_c("att_",
                                     as.character(.data$measured_att))) %>%
      dplyr::filter(.data$measured != 0) %>%
      dplyr::group_by(.data$profile, .data$item_id) %>%
      dplyr::mutate(meas =
                      stringr::str_c("att_",
                                     as.character(dplyr::row_number()))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"measured_att") %>%
      tidyr::pivot_wider(names_from = "meas", values_from = "measured") %>%
      dplyr::mutate(param =
                      dplyr::case_when(.data$master < 1 ~ NA_character_,
                                       .data$master == 1 ~
                                         as.character(glue::glue("l{item_id}_2",
                                                                 "{att_1}",
                                                                 "{att_2}"))
                                       )) %>%
      dplyr::select("profile", "item_id", "param")
  } else {
    profile_item_interactions <- tibble::tibble(profile =
                                                  rep(1:(2^ncol(q_matrix)),
                                                      each = nrow(q_matrix)),
                                                item_id =
                                                  rep(seq_len(nrow(q_matrix)),
                                                      times =
                                                      (2^ncol(q_matrix)))) %>%
      dplyr::mutate(param = NA_character_)
  }

  pi_mat <- calc_pi_mat(q_matrix, profs, profile_item_interactions, fungible)

  if (multithread) {
    stan_functions <- stan_func_block("TDCM", t, multithread, repeated_items)
    stan_functions <- stringr::str_replace_all(stan_functions, "\\{\\{", "{")
    stan_functions <- stringr::str_replace_all(stan_functions, "\\}\\}", "}")
  } else {
    stan_functions <- ""
  }

  stan_data <- data_block("TDCM", t, multithread)
  stan_data <- stringr::str_replace_all(stan_data, "\\{\\{", "{")
  stan_data <- stringr::str_replace_all(stan_data, "\\}\\}", "}")

  stan_transformed_data <- trans_data_block("TDCM", t)
  stan_transformed_data <- stringr::str_replace_all(stan_transformed_data,
                                                    "\\{\\{", "{")
  stan_transformed_data <- stringr::str_replace_all(stan_transformed_data,
                                                    "\\}\\}", "}")

  stan_parameters <- parameter_block("TDCM", t, int0, mef, int2)
  stan_parameters <- stringr::str_replace_all(stan_parameters, "\\{\\{", "{")
  stan_parameters <- stringr::str_replace_all(stan_parameters, "\\}\\}", "}")

  stan_transformed_parameters <- trans_param_block("TDCM", pi_mat, multithread)
  stan_transformed_parameters <-
    stringr::str_replace_all(stan_transformed_parameters, "\\{\\{", "{")
  stan_transformed_parameters <-
    stringr::str_replace_all(stan_transformed_parameters, "\\}\\}", "}")

  stan_model <- model_block("TDCM", t, priors, multithread, repeated_items)
  stan_model <- stringr::str_replace_all(stan_model, "\\{\\{", "{")
  stan_model <- stringr::str_replace_all(stan_model, "\\}\\}", "}")

  stan_generated_quantities <- gqs_block("TDCM", t, multithread, repeated_items)
  stan_generated_quantities <-
    stringr::str_replace_all(stan_generated_quantities, "\\{\\{", "{")
  stan_generated_quantities <-
    stringr::str_replace_all(stan_generated_quantities, "\\}\\}", "}")

  if (multithread) {
    stan_code <- list(functions = stan_functions,
                      data = stan_data,
                      transformed_data = stan_transformed_data,
                      parameters = stan_parameters,
                      transformed_parameters = stan_transformed_parameters,
                      model = stan_model,
                      generated_quantities = stan_generated_quantities)
  } else {
    stan_code <- list(data = stan_data,
                      parameters = stan_parameters,
                      transformed_parameters = stan_transformed_parameters,
                      model = stan_model,
                      generated_quantities = stan_generated_quantities)
  }


  return(stan_code)
}
