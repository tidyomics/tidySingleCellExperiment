data(pbmc_small)

test_that("auxiliary_colData() exposes reduced dimensions, column-aligned", {
    aux <- auxiliary_colData(pbmc_small)
    expect_s4_class(aux, "DataFrame")
    expect_identical(nrow(aux), ncol(pbmc_small))
    expect_identical(rownames(aux), colnames(pbmc_small))
    expect_true(all(c("PC_1", "tSNE_1") %in% colnames(aux)))
    expect_identical(
        unname(aux$PC_1),
        unname(reducedDim(pbmc_small, "PCA")[, 1]))
    expect_identical(
        unname(aux$tSNE_1),
        unname(reducedDim(pbmc_small, "TSNE")[, 1]))
})

test_that("auxiliary_colData() handles objects without reduced dimensions", {
    df <- pbmc_small
    reducedDims(df) <- list()
    aux <- auxiliary_colData(df)
    expect_identical(dim(aux), c(ncol(df), 0L))
    expect_identical(rownames(aux), colnames(df))
})

test_that("auxiliary_colData() falls back to indices without cell names", {
    df <- pbmc_small
    colnames(df) <- NULL
    aux <- auxiliary_colData(df)
    expect_identical(nrow(aux), ncol(df))
    expect_identical(rownames(aux), as.character(seq_len(ncol(df))))
})

test_that("auxiliary_colData() disambiguates shared dimension names", {
    df <- pbmc_small
    reducedDims(df) <- list(
        A=reducedDim(pbmc_small, "TSNE"),
        B=reducedDim(pbmc_small, "TSNE"))
    expect_false(anyDuplicated(colnames(auxiliary_colData(df))) > 0)
})

test_that("auxiliary_colData() honours n_dimensions_to_return", {
    expect_identical(
        ncol(auxiliary_colData(pbmc_small, n_dimensions_to_return=2)), 4L)
})

test_that("auxiliary names agree with the view-only columns", {
    # Dedup conventions differ for reductions sharing dimension names, so only
    # the unambiguous case is required to agree
    expect_setequal(
        get_special_columns(pbmc_small),
        colnames(auxiliary_colData(pbmc_small)))
    # mutate() still refuses to write a reduced dimension into colData
    expect_error(mutate(pbmc_small, PC_1=1), "view only")
})
