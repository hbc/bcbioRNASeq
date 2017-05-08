# bcbioRnaseq ====
library(bcbioRnaseq)
if (file.exists("data/run.rda")) {
    data(run)
} else {
    create_new_project()
    run <- load_run(
        file.path("upload_dir",
                  "YYYY-MM-DD_illumina_rnaseq",
                  "project-summary.yaml"),
        metadata_file = "meta/illumina_rnaseq.xlsx",
        intgroup = c("genotype", "treatment"))
    save_data(run)
}

# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))

# knitr ====
opts_chunk$set(
    autodep = TRUE,
    cache = TRUE,
    cache.lazy = TRUE,
    error = FALSE,
    fig.align = "center",
    fig.height = 8,
    fig.keep = "all",
    fig.path = "figures/",
    fig.width = 8,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)
