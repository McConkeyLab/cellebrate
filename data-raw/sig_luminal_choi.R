## code to prepare `sig_luminal_choi` dataset goes here
library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(homologene)
library(Seurat)
library(tibble)

choi_1 <- getGEO("GSE47993")[[1]]
choi_2 <- getGEO("GSE48124")[[1]]

get_rosi_samples <- function(x) {
  pd <- pData(x)
  get <- which(startsWith(pd$title, "UM-UC9") | startsWith(pd$title, "UC7"))
  x[, get]
}

rosi <- get_rosi_samples(combine(choi_1, choi_2))

get_rosi_degs <- function(x, coef, number, p.value, lfc) {
  mm <- model.matrix(~ c(rep(c("uc9", "uc7"), each = 6)) + c(rep(c("N", "Y", "Y", "N"), each = 3)))
  colnames(mm) <- c("intercept", "cell_line", "has_rosi")
  lmFit(x, mm) |>
    eBayes() |>
    topTable(coef = coef, number = number, p.value = p.value, lfc = lfc) |>
    as_tibble() |>
    select(
      symbol = Symbol,
      lfc = logFC,
      padj = adj.P.Val,
      mean_exp = AveExpr
    )
}

rosi_degs <- get_rosi_degs(rosi, coef = "has_rosi", number = Inf, p.value = 0.01, lfc = 1)

updated_names <- GeneSymbolThesarus(rosi_degs$symbol, search.types = "prev_symbol")

rosi_degs <- left_join(rosi_degs, enframe(updated_names), by = c("symbol" = "name"))

rosi_degs <- mutate(rosi_degs, symbol = ifelse(!is.na(value), value, symbol))

# These genes were found here:
# https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
# I chose only genes that had a 1-1 mapping
jax <- tibble(
  human = c("BHLE41", "CALHM3", "DHRS2", "CRIP2", "APOC1"),
  mouse = c("Bhle41", "Calhm3", "Dhrs2", "Crip2", "Apoc1")
)

converted_genes <- human2mouse(rosi_degs$symbol, db = homologeneData2)

gene_list <- tibble(
  human = converted_genes$humanGene,
  mouse = converted_genes$mouseGene
)

gene_list <- rbind(gene_list, jax)

# In case of duplicate, choose higher mean_exp
rosi_degs_no_dup <- rosi_degs |>
  group_by(symbol) |>
  mutate(keep = max(mean_exp) == mean_exp) |>
  filter(keep)

sig_luminal_choi <- gene_list |>
  filter(human != "HSD3B1") |> # 1 to many (H -> M) mapping
  right_join(rosi_degs_no_dup, by = c(human = "symbol")) |>
  select(-value, -keep)

usethis::use_data(sig_luminal_choi, overwrite = TRUE)
