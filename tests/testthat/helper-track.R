#' Generate a random track name for testing
#'
#' Creates a unique track name by combining "test." prefix with a random
#' alphanumeric string. Useful for avoiding name collisions in parallel tests.
#'
#' @param prefix Character string to use as prefix (default: "test")
#' @param length Integer, length of random suffix (default: 12)
#' @return Character string with format "prefix.randomstring"
#' @examples
#' random_track_name()
#' random_track_name(prefix = "mytrack", length = 8)
#' @noRd
random_track_name <- function(prefix = "temp", length = 12) {
    random_suffix <- paste(sample(c(letters, as.character(0:9)), length, replace = TRUE), collapse = "")
    paste0(prefix, ".", random_suffix)
}
