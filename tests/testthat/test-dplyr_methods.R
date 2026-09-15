data(pbmc_small)
df <- pbmc_small

# Only the verbs that tidySummarizedExperiment does not provide are defined in
# this package. Everything else is exercised by tidySummarizedExperiment's own
# tests; here we check that those verbs resolve to it and keep the object.

test_that("delegated verbs resolve to tidySummarizedExperiment", {
    delegated <- c("mutate", "filter", "select", "rename", "slice", "pull",
        "left_join", "inner_join", "right_join", "full_join", "distinct",
        "group_by", "summarise", "count", "rowwise", "sample_n", "sample_frac",
        "group_split", "bind_rows", "append_samples", "nest", "unnest",
        "extract", "unite", "separate", "pivot_longer", "plot_ly")
    # This package must not define them; the behavioural tests below and in
    # tidySummarizedExperiment check that they still work on an object
    for (v in delegated) {
        expect_null(
            utils::getS3method(v, "SingleCellExperiment", optional=TRUE),
            info=v)
    }
})

test_that("delegated verbs keep the object and its reduced dimensions", {
    fd <- filter(df, groups == "g1")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(reducedDimNames(fd), reducedDimNames(df))
    expect_true(ncol(fd) < ncol(df))

    fd <- mutate(df, blah=1)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_true("blah" %in% colnames(colData(fd)))
    expect_identical(dim(fd), dim(df))
})

test_that("delegated verbs never write reduced dimensions into colData", {
    y <- data.frame(.sample=colnames(df), ann=seq_len(ncol(df)))
    for (fd in list(
        suppressMessages(left_join(df, y, by=".sample")),
        mutate(df, blah=1),
        filter(df, groups == "g1"),
        tidyr::separate(df, groups, c("p", "q"), sep=1))) {
        expect_false(any(
            c("PC_1", "tSNE_1") %in% colnames(colData(fd))),
            info=paste(dim(fd), collapse="x"))
        # nor turn them into assays
        expect_identical(assayNames(fd), assayNames(df))
        expect_identical(reducedDimNames(fd), reducedDimNames(df))
    }
})

test_that("reduced dimensions are usable in delegated verbs", {
    fd <- filter(df, PC_1 > 0)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(ncol(fd), sum(reducedDim(df, "PCA")[, 1] > 0))
    # and are not written into colData
    fd <- mutate(df, pc_pos=PC_1 > 0)
    expect_true("pc_pos" %in% colnames(colData(fd)))
    expect_false("PC_1" %in% colnames(colData(fd)))
})

test_that("delegated append_samples() and nest()/unnest() keep the subclass", {
    # Binding an object to itself duplicates the cell names, which is warned
    fd <- suppressWarnings(append_samples(pbmc_small, pbmc_small))
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(ncol(fd), 2L*ncol(df))
    expect_identical(reducedDimNames(fd), reducedDimNames(df))

    fd <- unnest(nest(df, data=-groups), data)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(dim(fd), dim(df))
})

test_that("bind_cols() stays cell-aligned", {
    # tidySummarizedExperiment's bind_cols() binds against the
    # feature-by-sample tibble, so this method is kept here
    fd <- bind_cols(df, data.frame(x=seq_len(ncol(df))))
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(dim(fd), dim(df))
    expect_identical(fd$x, seq_len(ncol(df)))
})

test_that("arrange()", {
    fd <- arrange(df, nFeature_RNA)
    expect_s3_class(fd, "tbl_df")
    expect_identical(nrow(fd), nrow(df)*ncol(df))
    expect_true(!is.unsorted(fd$nFeature_RNA))
})

test_that("add_count()", {
    fd <- add_count(df, groups)
    expect_s3_class(fd, "tbl_df")
    expect_identical(nrow(fd), nrow(df)*ncol(df))
    expect_true("n" %in% colnames(fd))
    # 'n' counts rows of the feature-by-sample tibble within each group
    expect_identical(
        sort(unique(fd$n)),
        sort(unname(as.integer(table(as_tibble(df)$groups)))))
})

test_that("anti_join()", {
    y <- tibble::tibble(groups=unique(df$groups)[1])
    fd <- anti_join(df, y, by="groups")
    expect_s3_class(fd, "tbl_df")
    expect_false(any(fd$groups %in% y$groups))
})

test_that("slice_head() and slice_tail()", {
    expect_identical(nrow(slice_head(df, n=3)), 3L)
    expect_identical(nrow(slice_tail(df, n=3)), 3L)
    expect_s3_class(slice_head(df, n=1), "tbl_df")
})

test_that("slice_min() and slice_max()", {
    fd <- slice_min(df, nFeature_RNA, n=1, with_ties=FALSE)
    expect_identical(nrow(fd), 1L)
    expect_identical(fd$nFeature_RNA, min(as_tibble(df)$nFeature_RNA))
    fd <- slice_max(df, nFeature_RNA, n=1, with_ties=FALSE)
    expect_identical(fd$nFeature_RNA, max(as_tibble(df)$nFeature_RNA))
})

test_that("slice_sample()", {
    expect_identical(nrow(slice_sample(df, n=5)), 5L)
    expect_identical(nrow(slice_sample(df, prop=0)), 0L)
    expect_s3_class(slice_sample(df, n=1), "tbl_df")
})
