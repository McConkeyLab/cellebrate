#' Plot the expression of one or two genes with an optional stratifier
#'
#' @param dds A `DESeqDataSet` object
#' @param genes A character vector of one or two genes.
#' @param stratifier A name of a column of the `dds` `colData` that will be used
#'   to distinguish groups - see the dispatched function documentation for more
#'   details (see details)
#' @param assay integer. Assay position in the supplied `dds` that corresponds
#'   to normalized counts.
#' @return a ggplot
#' @details `eda_n_gene` will dispatch either `eda_one_gene` or `eda_two_gene`
#'   depending on the length of the `genes` argument.
#' @export
eda_n_gene <- function(dds, genes, stratifier = NULL, assay = 2) {
  if (length(genes) == 2) {
    eda_two_gene(dds, genes[1], genes[2], stratifier, assay)
  } else if (length(genes) == 1) {
    eda_one_gene(dds, genes[1], stratifier, assay)
  } else {
    stop("No suitable method to plot this number of genes.")
  }
}


#' Plot the expression of two genes
#'
#' @param dds A `DESeqDataSet` object
#' @param gene_x character. An HGNC symbol, like "KRT14", to be plotted on the
#'   x-axis
#' @param gene_y character. An HGNC symbol, like "SRC", to be plotted on the
#'   y-axis
#' @param stratifier A name of a column of the `dds` `colData` that will be used
#'   to color the points of the scatterplot
#' @param assay integer. Assay position in the supplied `dds` that corresponds
#'   to normalized counts.
#' @importFrom ggplot2 aes
#' @return a ggplot
#' @export
eda_two_gene <- function(dds,
                         gene_x,
                         gene_y,
                         stratifier = NULL,
                         assay = 2) {
  if (is.null(stratifier) &
      "clade" %in% colnames(SummarizedExperiment::colData(dds))) {
    stratifier <- "clade"
  }

  gene_inds <- get_gene_index(dds, c(gene_x, gene_y))
  gene_expr <- get_gene_expression(dds, gene_inds, assay)

  dds |>
    SummarizedExperiment::colData() |>
    as.data.frame() |>
    cbind(gene_expr) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[gene_x]],
                                 y = .data[[gene_y]],
                                 color = .data[[stratifier]])) +
    ggplot2::geom_point() +
    ggplot2::labs(x = gene_x, y = gene_y) +
    ggrepel::geom_text_repel(aes(label = .data$cell)) +
    bladdr::theme_tufte()
}

#' Plot the expression of one gene across a stratifying argument
#'
#' @param dds A `DESeqDataSet` object
#' @param gene A character vector HGNC Symbol, like "KRT14"
#' @param stratifier A name of a column of the `dds` `colData` that will be used
#'   as the x-axis to stratify gene expression into groups
#' @param assay The index of the assay that contains the normalized counts.
#' @param label character. Column to label points with
#' @importFrom ggplot2 aes
#' @return A ggplot
#' @export
eda_one_gene <- function(dds,
                         gene,
                         stratifier = NULL,
                         assay = 2,
                         label = NULL) {
  if (is.null(stratifier) &
      "clade" %in% colnames(SummarizedExperiment::colData(dds))) {
    stratifier <- "clade"
  }

  if (is.null(label) &
      "cell" %in% colnames(SummarizedExperiment::colData(dds))) {
    label <- "cell"
  }

  gene_ind <- get_gene_index(dds, gene)
  expression <- get_gene_expression(dds, gene_ind, assay)

  dds |>
    SummarizedExperiment::colData() |>
    as.data.frame() |>
    cbind(expression) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[stratifier]],
                                 y = .data[[gene]],
                                 color = .data[[stratifier]])) +
    ggplot2::geom_point() +
    ggplot2::labs(y = gene, x = NULL) +
    ggrepel::geom_text_repel(aes(label = .data[[label]]), nudge_x = 0.5) +
    bladdr::theme_tufte() +
    ggplot2::theme(legend.position = "none")
}

#' Get indices of genes by name in a dds
#' @noRd
get_gene_index <- function(dds, genes) {
  inds <- which(rownames(dds) %in% genes)
  if (identical(inds, integer(0))) {
    stop("None of the supplied genes exist in this dataset.")
  } else if (length(genes) != length(inds)) {
    found <- rownames(dds)[inds]
    missing <- setdiff(genes, found)
    warning(
      "Dataset does not contain gene: ",
      paste(missing, collapse = ", "),
      "\n"
    )
  }
  inds
}

#' Get expression of gene by indices in a dds
#'
#' @param dds A DESeqDataSet
#' @param gene_indices integer vector that refer to the row number of genes to retrieve
#' @param assay integer, assay index that refers to the desired expression data slot
#' @return A dataframe where row names are sample names, and column names are
#'   gene names
#' @details Importantly, returns data as a `data.frame` without dropping
#'   dimensions.
get_gene_expression <- function(dds, gene_indices, assay) {
  SummarizedExperiment::assay(dds, assay)[gene_indices, , drop = FALSE] |>
    t()
}
