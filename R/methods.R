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
                            exclude_zeros=exclude_zeros,
                            ...)) %>%
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
                    all=all,
                    ...))
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
#'   aggregate_cells(c(groups, ident), assays="counts")
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom ttservice aggregate_cells
#' @importFrom SummarizedExperiment assays assays<- assayNames
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stringr str_remove str_remove_all
#' @importFrom dplyr full_join
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by
#' @importFrom dplyr pick
#' @importFrom dplyr group_keys
#' @importFrom dplyr group_rows
#' @importFrom dplyr group_split
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr everything
#' @importFrom tidyr unite
#' @importFrom tidyr separate
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom purrr set_names
#' @importFrom purrr list_transpose
#' @importFrom S4Vectors DataFrame
#' @importFrom methods as
#'
setMethod("aggregate_cells", "SingleCellExperiment", function(.data,
    .sample=NULL, slot="data", assays=NULL,
    aggregation_function=Matrix::rowSums,
    ...) {
    # Fix NOTEs
    feature <- NULL
    .feature <- NULL
    my_id_to_split_by___ <- NULL
    assay_name <- NULL
    grouping_factor <- NULL
    assay_type <- NULL
    .sample <- enquo(.sample)

    arg_list <- c(mget(ls(environment(), sorted=FALSE)), match.call(expand.dots=FALSE)$...)
    assays_to_use <- eval(arg_list$assays)
    if (is.null(assays_to_use)) assays_to_use <- tail(names(assays(.data)), n=1)

    # Get information on sample groups
    sample_groups <- .data |>
        as_tibble() |>
        group_by(pick({{ .sample }}))

    sample_group_idx <- sample_groups |>
        group_rows()

    sample_group_keys <- sample_groups |>
        group_keys()

    .sample_names <- colnames(sample_group_keys)

    grouping_factor_names <- sample_group_keys |>
        unite(col="grouping_factor", !!.sample, sep="___") |>
        pull(grouping_factor)

    # Split sce object by groups
    sce_split <- map(.x=seq_along(sample_group_idx), .f=\(.num) .data[, sample_group_idx[[.num]]]) |>
        purrr::set_names(grouping_factor_names)

    grouping_factor <-
        .data |>
        colData() |>
        as_tibble() |>
        select(!!.sample) |>
        suppressMessages() |>
        unite("my_id_to_split_by___", !!.sample, sep="___") |>
        pull(my_id_to_split_by___) |>
        as.factor()

    # Add count of aggregated cells
    list_count_cells <- table(grouping_factor) |>
        enframe(name="grouping_factor", value=".aggregated_cells") |>
        mutate(.aggregated_cells=as.integer(.aggregated_cells))

    # Subset features based on selected assays
    feature_df <- get_all_features(.data)
    selected_features <- feature_df[feature_df$assay_id %in% assays_to_use, ]
    selected_experiments_list <- split(
        x=selected_features, f=as.character(selected_features$exp_id))
    if ("RNA" %in% names(selected_experiments_list))
        selected_experiments_list <- selected_experiments_list[
            c("RNA", setdiff(names(selected_experiments_list), "RNA"))]

    # Aggregate cells based on selected features from any assay / experiment type. Output is a tibble.
    aggregate_assays_fun <- function(exp) {
        # Check where the assay data needs to be taken from (main RNA experiment or altExp)
        selected_exp <- unique(exp$exp_id)
        selected_assays <- exp |> distinct(assay_name, .keep_all=TRUE)
        if (selected_exp == "RNA") {
            aggregate_sce_fun <- function(sce) {
                aggregated_vals <- assays(sce)[selected_assays$assay_name] |>
                    as.list() |>
                    map(.f=\(.list) aggregation_function(.list))
                map(.x=seq_along(aggregated_vals), \(.num) enframe(
                    x=aggregated_vals[[.num]],
                    name=".feature",
                    value=selected_assays$assay_id[[.num]])) |>
                    suppressMessages(reduce(full_join))
            }
            aggregated_list <- lapply(sce_split, aggregate_sce_fun) |>
                purrr::list_transpose() |>
                map(.f=\(.list) .list |> bind_rows(.id="grouping_factor"))
            interim_res <- map(.x=seq_along(aggregated_list), .f=\(.num) aggregated_list[[.num]] |>
                    separate(col=grouping_factor, into=.sample_names, sep="___")) |>
                purrr::set_names(nm=selected_exp)
            map(.x=seq_along(interim_res), .f=\(.num) interim_res[[.num]] |>
                    mutate(assay_type=names(interim_res)[[.num]])) |>
                purrr::reduce(full_join) |>
                select(assay_type, everything())
        } else {
            # aggregate from altExp
            aggregate_sce_fun <- function(sce) {
                aggregated_vals <- assays(altExps(sce)[[selected_exp]])[selected_assays$assay_name] |>
                    as.list() |>
                    purrr::set_names(selected_assays$assay_id) |>
                    map(.f=\(.list) aggregation_function(.list))
                map(.x=seq_along(aggregated_vals), \(.num) enframe(
                    x=aggregated_vals[[.num]],
                    name=".feature",
                    value=selected_assays$assay_id[[.num]])) |>
                    suppressMessages(reduce(full_join))
            }
            aggregated_list <- lapply(sce_split, aggregate_sce_fun) |>
                purrr::list_transpose() |>
                map(.f=\(.list) .list |> bind_rows(.id="grouping_factor"))
            interim_res <- map(.x=seq_along(aggregated_list), .f=\(.num) aggregated_list[[.num]] |>
                    separate(col=grouping_factor, into=.sample_names, sep="___")) |>
                purrr::set_names(nm=selected_exp)
            map(.x=seq_along(interim_res), .f=\(.num) interim_res[[.num]] |>
                    mutate(assay_type=names(interim_res)[[.num]])) |>
                purrr::reduce(full_join) |>
                select(assay_type, everything())
        }
    }

    # Join tibbles from each assay / experiment type into a single tibble.
    se <- lapply(selected_experiments_list, aggregate_assays_fun) |>
        purrr::reduce(full_join) |>
        suppressMessages()

    # Sometimes feature names can be duplicated in multiple assays, e.g. CD4 in RNA and ADT.
    any_feat_duplicated <- se |>
        distinct(assay_type, .feature) |>
        pull(.feature) |>
        duplicated() |>
        any()

    if (any_feat_duplicated) {
        warning("tidySingleCellExperiment says: The selected assays have overlapping feature names. The feature names have been combined with the selected assay_type, to keep the rownames of the SingleCellExperiment unique. You can find the original feature names in the feature_original column of the rowData slot of your object.")
        orig_features <- se |>
            distinct(assay_type, .feature)
        dup_features <- orig_features |>
            filter(duplicated(.feature)) |>
            pull(.feature)
        se <- se |>
            mutate(.feature=case_when(
                .feature %in% dup_features ~ str_c(assay_type, .feature, sep=".."),
                .default=.feature))
    }
    # Turn tibble into SummarizedExperiment object
    se <- se |>
        as_SummarizedExperiment(
            .sample=.sample_names,
            .transcript=.feature,
            .abundance=setdiff(colnames(se), c("assay_type", .sample_names, ".feature")))
    # Manually force the assay_type data to live in the rowData slot if it does not already
    if (exists("assay_type", where=as.data.frame(colData(se)))) {
        rowData(se) <- rownames(se) |>
            enframe(name=NULL, value="rowname") |>
            mutate(assay_type=unique(colData(se)$assay_type)) |>
            tibble::column_to_rownames() |>
            as.data.frame() |>
            as(Class="DataFrame")
        colData(se)$assay_type <- NULL
    }
    # Add original feature name information to the rowData slot
    if (any_feat_duplicated) {
        rowData(se) <- rowData(se) |>
            as.data.frame() |>
            tibble::rownames_to_column() |>
            mutate(feature_original=rowname,
                feature_original=str_remove_all(string=feature_original, pattern=".+(?=\\.\\.)"),
                feature_original=str_remove_all(string=feature_original, pattern="^\\..")) |>
            tibble::column_to_rownames() |>
            as(Class="DataFrame")
    }
    return(se)
})
