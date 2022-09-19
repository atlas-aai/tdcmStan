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
