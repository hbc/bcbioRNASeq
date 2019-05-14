# Set seed for reproducibility
set.seed(1454944673L)

requireNamespace(package = "knitr", quietly = TRUE)
knitr::opts_chunk[["set"]](
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

# Set default ggplot2 theme.
requireNamespace(package = "ggplot2", quietly = TRUE)
requireNamespace(package = "acidplots", quietly = TRUE)
ggplot2::theme_set(
    acidplots::theme_paperwhite(
        base_size = 14L,
        legend_position = "right"
    )
)
