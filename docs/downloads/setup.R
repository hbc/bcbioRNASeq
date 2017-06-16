# knitr ====
library(knitr)
opts_chunk$set(
    autodep = TRUE,
    bootstrap.show.code = FALSE,
    cache = TRUE,
    cache.lazy = TRUE,
    dev = c("png", "pdf", "svg"),
    error = TRUE,
    fig.height = 7,
    fig.retina = 2,
    fig.width = 7,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)



# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))



# bcbioRnaseq ====
library(bcbioRnaseq)
# Upload directory: location of the bcbio-nextgen final output
upload_dir <- "final"

# Columns in the colData in metadata you want to use
groups_of_interest <- c("genotype", "treatment")

# Output directory paths
out_path <- getwd()
data_out <- file.path(out_path, "data")
count_out <- file.path(out_path, "results", "counts")
res_out <- file.path(out_path, "results", "de")

# bcbioRnaDataSet
if (file.exists(file.path(data_out, "bcb.rda"))) {
    load(file.path(data_out, "bcb.rda"))
} else {
    bcb <- load_run(
        upload_dir = file.path(upload_dir),
        groups_of_interest = groups_of_interest)
    save_data(bcb, dir = data_out)
}

# DESeqDataSet
if (file.exists(file.path(data_out, "dds.rda")))
    load(file.path(data_out, "dds.rda"))
if (file.exists(file.path(data_out, "ddsde.rda")))
    load(file.path(data_out, "ddsde.rda"))
