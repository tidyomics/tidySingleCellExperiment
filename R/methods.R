#' @importFrom methods getMethod
setMethod(
    f="show",
    signature="SingleCellExperiment",
    definition=function(object) {
        opt <- getOption("restore_SingleCellExperiment_show", default=FALSE)
        if (isTRUE(opt)) {
            f <- getMethod(
                f="show",
                signature="SummarizedExperiment",
                where=asNamespace(ns="SummarizedExperiment"))
            f(object=object)
        } else { print(object) }
    }
)

setClass("tidySingleCellExperiment", contains="SingleCellExperiment")

#' @name join_features
#' @rdname join_features
#' @inherit ttservice::join_features
#' @aliases join_features,SingleCellExperiment-method
#'
#' @return A `tidySingleCellExperiment` object
#'   containing information for the specified features.
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small %>% join_features(
#'   features=c("HLA-DRA", "LYZ"))
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' @importFrom magrittr "%>%"
#' @importFrom dplyr contains
#' @importFrom dplyr everything
#' @importFrom ttservice join_features
#' @importFrom stringr str_c
#' @importFrom stringr str_subset
#' @export
setMethod("join_features", "SingleCellExperiment", function(.data,
    features=NULL, all=FALSE, exclude_zeros=FALSE, shape="wide", ...) {
    # CRAN Note
    .cell <- NULL
    .feature <- NULL

    # Shape is long
    if (shape == "long") {
      
        # Suppress generic data frame creation message produced by left_join
        suppressMessages({
            .data <-
                .data %>%
                    left_join(
                        by=c_(.data)$name,
                        get_abundance_sc_long(
                            .data=.data,
                            features=features,
                            all=all,
                            exclude_zeros=exclude_zeros)) %>%
                    select(!!c_(.data)$symbol, .feature,
                        contains(".abundance"), everything())
        })
      
        # Provide data frame creation and abundance column message
        if (any(class(.data) == "tbl_df")) {
            
            abundance_columns <-
                .data %>%
                colnames() %>%
                stringr::str_subset('.abundance_')
            
            message(stringr::str_c("tidySingleCellExperiment says: join_features produces",
                " duplicate cell names to accomadate the long data format. For this reason, a data", 
                " frame is returned for independent data analysis. Assay feature abundance is", 
                " appended as ", 
                stringr::str_flatten_comma(abundance_columns, last = " and "), "."
            ))
        }
      
        .data
        
    # Shape if wide
    } else {
        .data  %>%
            left_join(
                by=c_(.data)$name,
                get_abundance_sc_wide(
                    .data=.data,
                    features=features,
                    all=all, ...))
    }
})

#' @title (DEPRECATED) tidy for `SingleCellExperiment`
#' @name tidy.SingleCellExperiment
#'
#' @param x A `SingleCellExperiment` object.
#' @param ... Additional arguments passed to `generics::tidy`. (Unused.)
#' @return A `tidySingleCellExperiment` object. (DEPRECATED - not needed anymore)
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @importFrom generics tidy
#' @importFrom lifecycle deprecate_warn
#' @export
tidy.SingleCellExperiment <- function(x, ...) {

    # DEPRECATE
    deprecate_warn(
        when="1.1.1",
        what="tidy()",
        details="tidySingleCellExperiment says: tidy() is not needed anymore.")

    return(x)
}

#' @name aggregate_cells
#' @rdname aggregate_cells
#' @inherit ttservice::aggregate_cells
#' @aliases aggregate_cells,SingleCellExperiment-method
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small_pseudo_bulk <- pbmc_small |>
#'   aggregate_cells(.by = c(groups, ident), assays = "counts")
#'
#' @usage aggregate_cells(.data, slot = "data", assays = NULL,
#'   aggregation_function = Matrix::rowSums, .by = NULL, .sample = NULL, ...)
#'
#' @param .by Grouping columns (tidyverse-style).
#' @param .sample \lifecycle{soft-deprecated} Use `.by` instead.
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @return A `SummarizedExperiment` object with aggregated assays,
#'   preserving `rowData` and aggregated `colData`.
#'
#' @importFrom rlang as_quosure
#' @importFrom lifecycle deprecate_soft
#' @importFrom Matrix rowSums
#' @importFrom ttservice aggregate_cells
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#'
setMethod("aggregate_cells", "SingleCellExperiment", function(.data,
    .sample=NULL, slot="data", assays=NULL,
    aggregation_function=Matrix::rowSums,
    .by=NULL,
    ...) {

    # Re-capture NSE arguments from the original call, since S4 dispatch
    # evaluates arguments before passing to the method body. Walk up the call
    # stack past the S4 `.local` wrapper to find the actual user call.
    cl <- sys.call(sys.parent(1L))
    if (is.null(cl) || !is.call(cl)) cl <- sys.call()
    cl_names <- names(cl)

    # Determine which grouping argument was supplied (unevaluated)
    if (!is.null(cl_names) && ".by" %in% cl_names) {
        grouping_quo <- rlang::as_quosure(cl[[".by"]], parent.frame())
    } else if (!is.null(cl_names) && ".sample" %in% cl_names) {
        lifecycle::deprecate_soft(
            "1.20.2",
            "aggregate_cells(.sample=)",
            "aggregate_cells(.by=)",
            id = "aggregate_cells-sample-named"
        )
        grouping_quo <- rlang::as_quosure(cl[[".sample"]], parent.frame())
    } else if (length(cl) >= 3L) {
        # Positional second argument — treat as the deprecated .sample path
        lifecycle::deprecate_soft(
            "1.20.2",
            "aggregate_cells(.sample=)",
            "aggregate_cells(.by=)",
            id = "aggregate_cells-sample-positional"
        )
        grouping_quo <- rlang::as_quosure(cl[[3L]], parent.frame())
    } else {
        stop("`.by` must be specified.", call.=FALSE)
    }

    grouping_cols <- quosure_column_names(grouping_quo)

    use_assays <- if (is.null(assays)) assayNames(.data) else assays
    if (!all(use_assays %in% assayNames(.data))) {
        stop("assays not found in object", call.=FALSE)
    }

    if (!identical(aggregation_function, Matrix::rowSums)) {
        stop(
            "Only Matrix::rowSums is supported as aggregation_function.",
            call.=FALSE
        )
    }

    ids <- colData(.data)[, grouping_cols, drop=FALSE]

    check_and_install_packages("scrapper")

    aggregated_assays <- list()
    aggregated_col_data <- NULL

    for (i in seq_along(use_assays)) {
        at <- use_assays[[i]]
        agg <- scrapper::aggregateAcrossCells.se(
            .data,
            factors=ids,
            assay.type=at,
            output.prefix="",
            counts.name=if (i == 1L) "ncells" else NULL
        )
        aggregated_assays[[at]] <- assay(agg, "sums")
        if (i == 1L) {
            aggregated_col_data <- colData(agg)
        }
    }

    if (anyDuplicated(colnames(aggregated_col_data))) {
        aggregated_col_data <- aggregated_col_data[, !duplicated(colnames(aggregated_col_data)), drop=FALSE]
    }
    if ("ncells" %in% colnames(aggregated_col_data)) {
        colnames(aggregated_col_data)[colnames(aggregated_col_data) == "ncells"] <- ".aggregated_cells"
    }

    SummarizedExperiment::SummarizedExperiment(
        assays=aggregated_assays,
        rowData=rowData(.data),
        colData=aggregated_col_data
    )
})
