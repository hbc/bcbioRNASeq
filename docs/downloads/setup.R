# bcbioRnaseq ====
library(bcbioRnaseq)
out_path <- "." # where to save all output
data_out <- file.path(out_path, "data")
count_out <- file.path(out_path, "results", "counts")
res_out <- file.path(out_path, "results", "de")

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
