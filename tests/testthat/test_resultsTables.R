context("resultsTables")

load(system.file(
    file.path("extdata", "res.rda"),
    package = "bcbioRNASeq"))

test_that("resultsTables", {
    resTbl <- resultsTables(
        res,
        lfc = 0.25,
        summary = FALSE,
        write = FALSE)
    expect_identical(
        class(resTbl),
        "list"
    )
    tibble <- c("tbl_df", "tbl", "data.frame")
    expect_identical(
        lapply(resTbl, class),
        list("contrast" = "character",
             "alpha" = "numeric",
             "lfc" = "numeric",
             "all" = tibble,
             "deg" = tibble,
             "degLFC" = tibble,
             "degLFCUp" = tibble,
             "degLFCDown" = tibble)
    )
})
