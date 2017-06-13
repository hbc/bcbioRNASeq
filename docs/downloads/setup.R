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
    fig.retina = 2,
    fig.width = 7,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)

# parameters =====
upload_dir <- "" # bcbio_path_final_folder
genotype <- "" # variable in metadata you want to use, example group
out_path <- "." # where to save all output
# bcbioRnaseq ====
data_out <- file.path(out_path, "data")
count_out <- file.path(out_path, "results", "counts")
res_out <- file.path(out_path, "results", "de")

if (file.exists(file.path(data_out, "bcb.rda"))) {
    load(file.path(data_out, "bcb.rda"))
} else {
    bcb <- load_run_S4(
        file.path(upload_dir),
        intgroup = c(genotype))
    save_data(bcb, dir=data_out)
}

if (file.exists(file.path(data_out, "dds.rda")))
    load(file.path(data_out, "dds.rda"))

if (file.exists(file.path(data_out, "ddsde.rda")))
    load(file.path(data_out, "ddsde.rda"))
