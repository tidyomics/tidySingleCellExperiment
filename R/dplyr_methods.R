# Most dplyr verbs are provided by tidySummarizedExperiment, which operates on
# the feature-by-sample tibble abstraction. Only the verbs that
# tidySummarizedExperiment does not implement are defined here. They convert to
# the tibble abstraction and return a data frame, as the other verbs that cannot
# map their result back onto the object do.

#' @name arrange
#' @rdname arrange
#' @inherit dplyr::arrange
#' @family single table verbs
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> 
#'     arrange(nFeature_RNA)
#'     
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#'     
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange
#' @export
arrange.SingleCellExperiment <- function(.data, ..., .by_group=FALSE) {
    tidy_message(data_frame_returned_message)

    .data |>
        as_tibble() |>
        dplyr::arrange(..., .by_group=.by_group)
}

#' @name anti_join
#' @rdname anti_join
#' @inherit dplyr::anti_join
#'
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' tt |> anti_join(tt |> distinct(groups) |> mutate(new_column=1) |> slice(1))
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#'
#' @importFrom dplyr anti_join
#' @export
anti_join.SingleCellExperiment <- function(x, y, by=NULL, copy=FALSE, ...) {
    tidy_message(data_frame_returned_message)

    x |>
        as_tibble() |>
        dplyr::anti_join(y, by=by, copy=copy, ...)
}

#' @name add_count
#' @rdname add_count
#' @inherit dplyr::add_count
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> add_count(groups)
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#'
#' @importFrom dplyr add_count
#' @importFrom rlang enquo
#' @export
add_count.SingleCellExperiment <- function(x, ..., wt=NULL, sort=FALSE,
    name=NULL) {
    tidy_message(data_frame_returned_message)

    x |>
        as_tibble() |>
        dplyr::add_count(..., wt=!!enquo(wt), sort=sort, name=name)
}

#' @name slice_head
#' @rdname slice_head
#' @inherit dplyr::slice_head
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> slice_head(n=1)
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#'
#' @importFrom dplyr slice_head
#' @export
slice_head.SingleCellExperiment <- function(.data, ..., n, prop, by=NULL) {
    tidy_message(data_frame_returned_message)

    .data |>
        as_tibble() |>
        dplyr::slice_head(..., n=n, prop=prop, by={{ by }})
}

#' @rdname slice_head
#' @importFrom dplyr slice_tail
#' @export
slice_tail.SingleCellExperiment <- function(.data, ..., n, prop, by=NULL) {
    tidy_message(data_frame_returned_message)

    .data |>
        as_tibble() |>
        dplyr::slice_tail(..., n=n, prop=prop, by={{ by }})
}

#' @rdname slice_head
#' @importFrom dplyr slice_min
#' @export
slice_min.SingleCellExperiment <- function(.data, order_by, ..., n, prop,
    by=NULL, with_ties=TRUE, na_rm=FALSE) {
    tidy_message(data_frame_returned_message)

    .data |>
        as_tibble() |>
        dplyr::slice_min({{ order_by }}, ..., n=n, prop=prop, by={{ by }},
            with_ties=with_ties, na_rm=na_rm)
}

#' @rdname slice_head
#' @importFrom dplyr slice_max
#' @export
slice_max.SingleCellExperiment <- function(.data, order_by, ..., n, prop,
    by=NULL, with_ties=TRUE, na_rm=FALSE) {
    tidy_message(data_frame_returned_message)

    .data |>
        as_tibble() |>
        dplyr::slice_max({{ order_by }}, ..., n=n, prop=prop, by={{ by }},
            with_ties=with_ties, na_rm=na_rm)
}

#' @rdname slice_head
#' @importFrom dplyr slice_sample
#' @export
slice_sample.SingleCellExperiment <- function(.data, ..., n, prop, by=NULL,
    weight_by=NULL, replace=FALSE) {
    tidy_message(data_frame_returned_message)

    # 'n' and 'prop' are left missing rather than NULL, as dplyr rejects
    # receiving both
    .data |>
        as_tibble() |>
        dplyr::slice_sample(..., n=n, prop=prop, by={{ by }},
            weight_by={{ weight_by }}, replace=replace)
}
# 'bind_rows()' and 'append_samples()' are provided by
# tidySummarizedExperiment, which column-binds the objects and carries the
# same deprecation of 'bind_rows()' in favour of 'append_samples()'.

#' @name bind_cols
#' @rdname bind_cols
#' @inherit ttservice::bind_cols
#'
#' @description
#' Appends cell-aligned columns to \code{colData()}. This method is kept in
#' this package because tidySummarizedExperiment's \code{bind_cols()} binds
#' against the feature-by-sample tibble, whose row count is a multiple of the
#' number of cells.
#'
#' @examples
#' data(pbmc_small)
#' # 'ttservice::bind_cols()' is the generic; dplyr's is not one
#' ttservice::bind_cols(pbmc_small, data.frame(x=seq_len(ncol(pbmc_small))))
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166-1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom rlang dots_values
#' @importFrom ttservice bind_cols
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
bind_cols_ <- function(..., .id=NULL) {
    tts <- tts <- flatten_if(dots_values(...), is_spliced)
    
    colData(tts[[1]]) <- bind_cols(colData(tts[[1]]) %>% as.data.frame(),
        tts[[2]], .id=.id) %>% DataFrame()
    
    tts[[1]]
}

#' @rdname bind_cols
#' @export
bind_cols.SingleCellExperiment <- bind_cols_

