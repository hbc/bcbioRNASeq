context("resultsTables")

test_that("resultsTables", {
    resTbl <- resultsTables(res, lfc = 0.25, write = FALSE)
    expect_equal(
        class(resTbl),
        "list"
    )
    expect_equal(
        names(resTbl),
        c("contrast",
          "alpha",
          "lfc",
          "all",
          "deg",
          "degLFC",
          "degLFCUp",
          "degLFCDown",
          "allFile",
          "degFile",
          "degLFCUpFile",
          "degLFCDownFile")
    )
})
