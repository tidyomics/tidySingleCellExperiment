.onAttach <- function(libname, pkgname) {
    attached <- tidyverse_attach()
    
    tidyprint::tidy_message(
        "By default SingleCellExperiment uses the standard display. For a tidy tibble-style display, run tidy_print_on(remember = TRUE).",
        type = "info"
    )
}
