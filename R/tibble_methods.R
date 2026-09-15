# 'as_tibble()' is provided by tidySummarizedExperiment. The tibble is
# feature-by-sample, keyed by '.feature' and '.sample', with reduced
# dimensions contributed as view-only columns via 'auxiliary_colData()'.

#' @name glimpse
#' @rdname glimpse
#' @inherit pillar::glimpse
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> glimpse()
#' 
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#' 
#' @importFrom tibble glimpse
#' @export
glimpse.tidySingleCellExperiment <- function(x, width=NULL, ...){
    x %>%
        as_tibble() %>%
        tibble::glimpse(width=width, ...)
}
