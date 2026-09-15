#' Objects exported from other packages
#'
#' The `auxiliary_colData()` generic is re-exported from
#' \pkg{tidySummarizedExperiment} so that it is available when only
#' \pkg{tidySingleCellExperiment} is attached. The import binding is kept, so
#' the method below registers against the generic's own S3 table.
#'
#' @name reexports
#' @keywords internal
#'
#' @importFrom tidySummarizedExperiment auxiliary_colData
#' @export
tidySummarizedExperiment::auxiliary_colData

#' Reduced dimensions as column-aligned auxiliary metadata
#'
#' Exposes `reducedDims()` to the tidy verbs provided by
#' \pkg{tidySummarizedExperiment} as view-only, column-aligned columns.
#' Each reduced dimension matrix contributes one column per dimension, so
#' e.g. `filter(sce, PC_1 > 0)` can be resolved without materialising the
#' feature-by-sample tibble. These columns are never written back to
#' `colData()`.
#'
#' @param x A `SingleCellExperiment`.
#' @param ... Passed to `get_special_datasets()`, e.g.
#'   `n_dimensions_to_return` to cap the number of dimensions per reduction.
#'
#' @return A [S4Vectors::DataFrame] with one row per column of `x`.
#'
#' @examples
#' data(pbmc_small)
#' auxiliary_colData(pbmc_small)[, seq_len(2)]
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @importFrom S4Vectors DataFrame
#' @export
auxiliary_colData.SingleCellExperiment <- function(x, ...) {

    cell_names <- colnames(x)
    if (is.null(cell_names)) cell_names <- as.character(seq_len(ncol(x)))

    mats <- get_special_datasets(x, ...)
    if (!length(mats)) return(DataFrame(row.names=cell_names))

    df <- do.call(cbind, lapply(mats, as.data.frame))

    # 'cbind()' prefixes with the reduction name; restore the dimension names
    # and disambiguate reductions that share dimension names
    colnames(df) <- make.unique(
        unlist(lapply(mats, colnames), use.names=FALSE))
    rownames(df) <- cell_names

    DataFrame(df, check.names=FALSE)
}
