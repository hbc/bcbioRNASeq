# bcbioRnaseq ====
library(bcbioRnaseq)
output_dir <- getwd()  # where to save all output
data_dir <- file.path(output_dir, "data")
count_dir <- file.path(output_dir, "results", "counts")
res_dir <- file.path(output_dir, "results", "de")

# knitr ====
library(knitr)
opts_chunk$set(
    autodep = TRUE,
    bootstrap.show.code = FALSE,
    cache = TRUE,
    cache.lazy = TRUE,
    dev = c("png", "pdf", "svg"),
    error = TRUE,
    fig.height = 6,
    fig.retina = 2,
    fig.width = 6,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)

# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))
