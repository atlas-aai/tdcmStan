check_multi_attribute <- function(x) {
  num_multi_att_items <- x %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(dplyr::c_across(dplyr::where(is.numeric)))) %>%
    dplyr::filter(.data$total > 1) %>%
    nrow()

  if (num_multi_att_items > 0) {
    stop(paste("Q-matrix cannot have items measuring multiple attributes with",
               "a fungible model."),
         call. = FALSE)
  } else {
    x
  }
}
