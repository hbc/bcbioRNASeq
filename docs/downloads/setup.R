library(knitr)
library(ggplot2)
library(bcbioRNASeq)

opts_chunk$set(
    autodep = TRUE,
    bootstrap.show.code = FALSE,
    cache = TRUE,
    cache.lazy = TRUE,
    dev = c("png", "pdf"),
    error = TRUE,
    fig.height = 8,
    fig.retina = 2,
    fig.width = 10,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)

theme_set(theme_light(base_size = 14))
theme_update(legend.position = "bottom")
