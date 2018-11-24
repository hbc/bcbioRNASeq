library(knitr)
library(basejump)
library(bcbioRNASeq)
library(DESeq2)
library(tidyverse)

# Set seed for reproducibility.
set.seed(1454944673L)

opts_chunk[["set"]](
    autodep = TRUE,
    bootstrap.show.code = FALSE,
    cache = TRUE,
    cache.lazy = TRUE,
    dev = c("png", "pdf"),
    fig.height = 10L,
    fig.retina = 2L,
    fig.width = 10L,
    highlight = TRUE,
    prompt = TRUE,
    tidy = FALSE
)

theme_set(
    theme_paperwhite(
        base_size = 14L,
        legend_position = "right"
    )
)
