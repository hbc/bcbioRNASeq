context("plotGene")

symbol <- rowData(bcb)[["symbol"]][1:3]
ensgene <- rowData(bcb)[["ensgene"]][1:3]

# Check for Rgs20, Rab23, Ncoa2 inside the ggplot objects

test_that("Gene symbol", {
    p <- plotGene(bcb, gene = symbol, format = "symbol")
    expect_is(p, "ggplot")
})

test_that("Ensembl gene identifier", {
    p <- plotGene(bcb, gene = ensgene, format = "ensgene")
    expect_is(p, "ggplot")
})

test_that("plotlist return", {
    list <- plotGene(
        bcb,
        gene = symbol,
        format = "symbol",
        returnList = TRUE)
    expect_is(list, "list")
    expect_true(
        lapply(list, function(x) is(x, "ggplot")) %>%
            unlist() %>%
            all()
    )
})
