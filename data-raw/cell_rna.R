## code to prepare `cell_rna` dataset goes here
library(bladdr)

cell_rna <- get_gbci_file(
  "Datasets/Cell Line RNA Seq Data/thirty-cell-lines-rnaseq/cell-lines_norm_clades.Rds"
) |>
    readRDS()

usethis::use_data(cell_rna, overwrite = TRUE)
