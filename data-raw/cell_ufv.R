## code to prepare `cell_dna` dataset goes here
library(httr2)
library(tidyverse)
library(fs)

cell_names <- c("uc5", "tccsup", "rt4", "uc12", "uc14", "uc6", "j82", "5637",
                "uc3", "rt112", "scaber", "control")

analysis_names <- paste0("2023-01-25_", cell_names)

# Should take about 10 seconds.
all_analyses <- request("https://ionreporter.thermofisher.com/api/v1/analysis") |>
    req_headers(`Content-Type` = "application/x-www-form-urlencoded",
                Authorization = Sys.getenv("IONREPORTER_KEY", )) |>
    req_url_query(format = "json",
                  view = "summary") |>
    req_perform() |>
    resp_body_json() |>
    lapply(\(x) enframe(x) |> pivot_wider()) |>
    bind_rows() |>
    unnest(cols = everything(), keep_empty = TRUE) |>
    unnest(cols = everything(), keep_empty = TRUE) # Needed twice for deeper lists

cell_analyses <- all_analyses |>
  filter(str_detect(samples, paste0(cell_names, collapse = "|")))

# Takes around 5 minutes or so.
get_variants_links <- function(analysis_name) {
  response <- request("https://ionreporter.thermofisher.com/api/v1/analysis") |>
    req_headers(`Content-Type` = "application/x-www-form-urlencoded",
                Authorization = Sys.getenv("IONREPORTER_KEY", )) |>
    req_url_query(format = "json",
                  type = "analysis",
                  name = analysis_name,
                  exclude = "reports") |>
    req_perform() |>
    resp_body_json()
  response[[1]]$data_links
}

linked_cell_analyses <- cell_analyses |>
  rowwise() |>
  mutate(links = list(get_variants_links(name)))

download_file <- function(url, destfile) {
  download.file(url, destfile, method = "libcurl")
  destfile
}

ufv_locs <- download_file(
  linked_cell_analyses$links |> lapply(`[[`, "unfiltered_variants") |> unlist(),
  path(tempdir(),
       linked_cell_analyses$name,
       ext = "zip"))

extract_variant_file <- function(zip_path, og_sample, variants_dir) {
  # IonReporter automatically does this w sample names
  og_sample <- str_replace_all(og_sample, " ", "_")
  unzip(zip_path,
        files = path("Variants", og_sample,
                         paste0(og_sample, "-full"), ext = "tsv"),
        exdir = variants_dir,
        junkpaths = TRUE)
  path(variants_dir, paste0(og_sample, "-full"), ext = "tsv")
}

ufv_extracted_locs <- mapply(extract_variant_file, ufv_locs, linked_cell_analyses$samples,
                             tempdir())

read_variant_call_file <- function(path, workflow, sample_name, og_sample) {
  data <- read_tsv(path, skip = 2, col_types = "ccciccddciccccc", show_col_types = FALSE)
  data$workflow <- workflow
  data$sample <- sample_name
  data$og_sample <- og_sample
  data <- relocate(data, sample)
  data
}

ufv <- mapply(read_variant_call_file,
              ufv_extracted_locs,
              "gbci62",
              linked_cell_analyses$name,
              linked_cell_analyses$samples,
              SIMPLIFY = FALSE) |>
  bind_rows()


get_unique_genes <- function(genes) {
  split_genes <- strsplit(genes, "\\|")[[1]]
  rm_pipes <- split_genes[split_genes != "|"]
  unique_genes <- unique(rm_pipes)
}

tidy_variants <- function(filtered_variants) {
  filtered_variants |>
    dplyr::select(-cnv_pvalue) |>
    rowwise() |>
    mutate(gene = list(get_unique_genes(gene)),
           first_gene = gene[1]) |>
    rename(cell_line = og_sample,
           locus = `# locus`,
           gene_aliases = gene,
           gene = first_gene,
           allele_freq_pct = `allele_frequency_%`,
           gene_function = `function`) |>
    mutate(cell_line = str_remove(cell_line, "_v1$")) |>
    select(-sample) |>
    relocate(cell_line, gene, type, pvalue, phred_qual_score)
}

cell_ufv <- tidy_variants(ufv)

usethis::use_data(cell_ufv, overwrite = TRUE)
