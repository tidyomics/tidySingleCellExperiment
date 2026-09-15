\name{NEWS}
\title{News for Package \pkg{tidySingleCellExperiment}}

\section{Changes in version 2.0.0}{
\itemize{
    \item Printing is now handled by \pkg{tidyprint}, matching tidySummarizedExperiment. The standard SingleCellExperiment display is the default; use \code{tidy_print_on()}/\code{tidy_print_off()} to toggle the tidy tibble-style display. The custom SingleCellExperiment printer has been removed.
    \item Package messages now use \code{tidyprint::tidy_message()} for consistent tidyomics styling.
    \item Added an \code{auxiliary_colData()} method for SingleCellExperiment, implementing the extension point from \pkg{tidySummarizedExperiment}. Reduced dimensions are now exposed as view-only, column-aligned columns, so tidySummarizedExperiment can resolve expressions such as \code{PC_1 > 0} against them without writing them back to \code{colData()}.
    \item Breaking: the tibble abstraction is now the one provided by \pkg{tidySummarizedExperiment}. \code{as_tibble()} returns one row per feature and cell, keyed by \code{.feature} and \code{.sample}, with assays as columns, rather than one row per cell keyed by \code{.cell}. Cells are the columns of the object and are therefore keyed by \code{.sample}; \code{.cell} no longer exists. This makes the abstraction consistent with the printed output and with the rest of tidyomics. Code that refers to \code{.cell} must be updated to \code{.sample}.
    \item Reduced dimensions now appear as columns of \code{as_tibble()}, contributed through \code{auxiliary_colData()}. They remain view-only and are never written back to \code{colData()}.
    \item Assay values are now available to the tidy verbs, so expressions such as \code{filter(counts > 0)} and \code{mutate(logged=log1p(counts))} work. This was not possible with the cell-wise abstraction.
    \item Removed the \pkg{dplyr}, \pkg{tidyr}, \pkg{tibble} and \pkg{plotly} methods that \pkg{tidySummarizedExperiment} already provides: \code{mutate()}, \code{filter()}, \code{select()}, \code{rename()}, \code{slice()}, \code{pull()}, the four mutating joins, \code{distinct()}, \code{group_by()}, \code{summarise()}, \code{summarize()}, \code{count()}, \code{rowwise()}, \code{sample_n()}, \code{sample_frac()}, \code{group_split()}, \code{nest()}, \code{unnest()}, \code{extract()}, \code{unite()}, \code{separate()}, \code{pivot_longer()}, \code{plot_ly()} and \code{ggplot()}. They now resolve to \pkg{tidySummarizedExperiment} and operate on the shared abstraction. This removes roughly 1,700 lines of duplicated implementation.
    \item Also removed \code{bind_rows()} and \code{append_samples()}, which \pkg{tidySummarizedExperiment} implements identically, including the deprecation of \code{bind_rows()} in favour of \code{append_samples()}.
    \item \code{bind_cols()} is kept, as it appends cell-aligned columns to \code{colData()}, whereas \pkg{tidySummarizedExperiment}'s binds against the feature-by-sample tibble.
    \item \code{arrange()}, \code{anti_join()}, \code{add_count()} and the \code{slice_*()} family are not provided by \pkg{tidySummarizedExperiment} and are kept here. They now return a data frame, as the other verbs that cannot map their result back onto the object do.
    \item \code{join_features(shape="long")} is keyed by \code{.sample} rather than \code{.cell}; the returned columns are otherwise unchanged.
    \item \code{join_features(shape="wide")} now keeps the original feature names, for example \code{HLA-DRB1}, where previously they were made syntactically valid, as in \code{HLA.DRB1}. Refer to such columns with backticks or \code{.data[["HLA-DRB1"]]}.
    \item Removed the deprecated \code{join_transcripts()}; use \code{join_features()} instead.
    \item \code{\%>\%} is no longer re-exported. It remains available because the package attaches \pkg{dplyr}, or load \pkg{magrittr} directly. The native pipe \code{|>} is recommended.
}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item Soft deprecated \code{bind_rows()} in favor of \code{append_samples()} from ttservice.
    \item Added \code{append_samples()} method for SingleCellExperiment objects.
    \item \code{bind_rows()} is not a generic method in dplyr and may cause conflicts.
    \item Users are encouraged to use \code{append_samples()} instead.
}}

\section{Changes in version 1.4.0, Bioconductor 3.14 Release}{
\itemize{
    \item Improved sample_n, and sample_frac functions.
    \item Add join_features prefix.
    \item Dropped tidy method as never needed.
    \item Add unnest_tidySingleCellExperiment for nested data that was not produce with tidySingleCellExperiment::nest() but rather with tidyr::nest().
}}

\section{Changes in version 1.5.1, Bioconductor 3.15 Release}{
\itemize{
    \item Rely of ttservice package for shared function with tidySingleCellExperiment to avoid clash
    \item Use .cell for cell column name to avoid errors when cell column is defined by the user
}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item \strong{BREAKING CHANGE}: Changed default shape parameter in \code{join_features()} from "long" to "wide". 
    This means that \code{join_features()} now returns a SingleCellExperiment object by default instead of a tibble. 
    To get the old behavior, explicitly specify \code{shape="long"}.
}}

