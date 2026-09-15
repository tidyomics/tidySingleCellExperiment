data(pbmc_small)
df <- pbmc_small

test_that("show() uses standard display by default", {
    options(tidyprint.use_tidy_print = FALSE)
    txt <- capture.output(show(df))
    expect_true(any(grepl("^class: SingleCellExperiment", txt)))
    expect_true(any(grepl("^dim: ", txt)))
})

test_that("show() uses tidyprint when enabled", {
    options(tidyprint.use_tidy_print = TRUE)
    on.exit(options(tidyprint.use_tidy_print = FALSE), add=TRUE)
    txt <- capture.output(show(df))
    expect_true(any(grepl("SummarizedExperiment-tibble abstraction", txt, fixed=TRUE)))
})

test_that("join_features() wide appends abundances per cell", {
    gs <- sample(rownames(df), 3)
    fd <- join_features(df, gs, assay="counts")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_null(fd$.feature)
    # The requested features become cell-aligned columns in colData, keeping
    # their original names
    expect_true(all(gs %in% colnames(colData(fd))))
    # Index the DataFrame directly, as as.data.frame() would mangle names
    # such as 'HLA-DRB1'
    expect_identical(
        unname(t(as.matrix(colData(fd)[, gs, drop=FALSE]))),
        as.matrix(unname(counts(df)[gs, ])))
    # Reduced dimensions stay view-only: neither copied into colData nor
    # turned into assays
    expect_identical(
        setdiff(colnames(colData(fd)), colnames(colData(df))), gs)
    expect_identical(assayNames(fd), assayNames(df))
    expect_identical(reducedDimNames(fd), reducedDimNames(df))
})

test_that("join_features() long returns one row per feature and cell", {
    gs <- sample(rownames(df), 3)
    fd <- join_features(df, gs, shape="long")
    expect_s3_class(fd, "tbl_df")
    expect_setequal(unique(fd$.feature), gs)
    expect_identical(nrow(fd), length(gs)*ncol(df))
    expect_true(all(table(fd$.feature) == ncol(df)))
    # Abundances match the assay, matched on feature and cell
    expected <- vapply(seq_len(nrow(fd)),
        \(i) counts(df)[fd$.feature[i], fd$.sample[i]],
        numeric(1))
    expect_equal(fd$.abundance_counts, expected)
})

test_that("as_tibble() is the feature-by-sample abstraction", {
    fd <- as_tibble(df)
    expect_s3_class(fd, "tbl_df")
    # One row per feature and sample (cell), not one row per cell
    expect_equal(nrow(fd), nrow(df)*ncol(df))
    expect_true(all(c(".feature", ".sample") %in% colnames(fd)))
    expect_false(".cell" %in% colnames(fd))
    # Assays are columns
    expect_true(all(assayNames(df) %in% colnames(fd)))
})

test_that("as_tibble() exposes reduced dimensions as view-only columns", {
    fd <- as_tibble(df)
    expect_true(all(c("PC_1", "tSNE_1") %in% colnames(fd)))
    # Column-aligned: constant within a cell
    n <- dplyr::n_distinct(fd[fd$.sample == fd$.sample[1], ][["PC_1"]])
    expect_identical(n, 1L)
    # Still absent from colData
    expect_false("PC_1" %in% colnames(colData(df)))
})

test_that("aggregate_cells()", {
    df$factor <- sample(gl(3, 1, ncol(df)))
    df$string <- sample(c("a", "b"), ncol(df), TRUE)
    tbl <- distinct(select(df, factor, string))
    fd <- aggregate_cells(df, c(factor, string))
    expect_identical(assayNames(fd), assayNames(df))
    # [HLC: aggregate_cells() currently
    # reorders features alphabetically]
    fd <- fd[rownames(df), ]
    expect_s4_class(fd, "SummarizedExperiment")
    expect_equal(dim(fd), c(nrow(df), nrow(tbl)))
    foo <- mapply(
        f=tbl$factor,
        s=tbl$string,
        \(f, s) {
            expect_identical(
                df |> 
                    filter(factor == f, string == s) |>
                    assay() |> rowSums() |> as.vector(),
                fd[, fd$factor == f & fd$string == s] |>
                    assay() |> as.vector())
        })
    # specified 'assays' are subsetted
    expect_error(aggregate_cells(df, c(factor, string), assays="x"))
    fd <- aggregate_cells(df, c(factor, string), assays="counts")
    expect_identical(assayNames(fd), "counts")
})
