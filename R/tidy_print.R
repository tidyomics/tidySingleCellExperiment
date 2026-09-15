#' Enable or disable tidy print for SingleCellExperiment
#'
#' Re-exports from \pkg{tidyprint} for convenience. When tidy print is enabled,
#' SingleCellExperiment objects use \pkg{tidyprint}'s compact tibble-style display.
#'
#' @importFrom tidyprint tidy_print_on tidy_print_off tidy_print_enabled
#' @importFrom tidyprint tidy_message
#'
#' @param remember If \code{TRUE}, saves the setting to a cache file so it
#'   persists across R sessions. If \code{FALSE} (default), the setting only
#'   affects the current session.
#'
#' @return \code{tidy_print_on()} and \code{tidy_print_off()} return \code{TRUE}
#'   or \code{FALSE} invisibly. \code{tidy_print_enabled()} returns a logical.
#'
#' @seealso \code{\link[tidyprint:tidy_print_toggle]{tidyprint::tidy_print_on}}
#'
#' @name tidy_print
#' @export
tidy_print_on <- tidyprint::tidy_print_on

#' @rdname tidy_print
#' @export
tidy_print_off <- tidyprint::tidy_print_off

#' @rdname tidy_print
#' @export
tidy_print_enabled <- tidyprint::tidy_print_enabled

#' Emit a warning with tidyprint-style display
#' @keywords internal
#' @noRd
tidy_warning <- function(message) {
    tidyprint::tidy_message(message, type = "warning")
    warning(message, call. = FALSE)
}

#' Emit an error with tidyprint-style display
#' @keywords internal
#' @noRd
tidy_stop <- function(message) {
    tidyprint::tidy_message(message, type = "danger")
    stop(message, call. = FALSE)
}
