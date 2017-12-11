context("plotGene")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

ensgene <- rowData(bcb)[["ensgene"]][1:3]
symbol <- rowData(bcb)[["symbol"]][1:3]

test_that("Ensembl gene identifier", {
    p <- plotGene(bcb, genes = ensgene, format = "ensgene")
    expect_is(p, "ggplot")
})

test_that("Gene symbol", {
    p <- plotGene(bcb, genes = symbol, format = "symbol")
    expect_is(p, "ggplot")
})

test_that("plotlist", {
    list <- plotGene(
        bcb,
        genes = ensgene,
        format = "ensgene",
        returnList = TRUE)
    expect_is(list, "list")
    expect_true(
        lapply(list, function(x) is(x, "ggplot")) %>%
            unlist() %>%
            all()
    )
})
