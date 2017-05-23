# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))

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
    fig.keep = "all",
    fig.path = "figures/",
    fig.retina = 2,
    fig.width = 7,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)

# bcbioRnaseq ====
library(bcbioRnaseq)
if (file.exists("data/run.rda")) {
    data(run)
} else {
    run <- load_run(
        file.path("upload_dir"),
        intgroup = c("genotype"))
    save_data(run)
}
