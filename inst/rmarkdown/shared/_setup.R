# Set seed for reproducibility.
set.seed(1454944673L)

knitr::opts_chunk[["set"]](
    autodep = TRUE,
    bootstrap.show.code = FALSE,
    # Enable caching with caution.
    cache = FALSE,
    cache.lazy = FALSE,
    comment = "",
    dev = c("png", "pdf"),
    fig.height = 10L,
    fig.retina = 2L,
    fig.width = 10L,
    highlight = TRUE,
    prompt = TRUE,
    tidy = FALSE
)

ggplot2::theme_set(
    basejump::theme_paperwhite(
        base_size = 14L,
        legend_position = "right"
    )
)
