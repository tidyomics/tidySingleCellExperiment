# 'ggplot()' is provided by tidySummarizedExperiment, which builds the plot
# from the 'as_tibble()' generic and therefore receives the cell-wise
# abstraction for SingleCellExperiment.

# addressing R CMD CHECK NOTE "no visible global function definition for 'aes'"
globalVariables("aes", "tidySingleCellExperiment")
