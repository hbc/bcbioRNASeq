context("plotGene")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

genes <- rownames(bcb)[1L:3L]

test_that("Ensembl gene identifiers", {
    p <- plotGene(bcb, genes = genes)
    expect_is(p, "ggplot")
})

test_that("plotlist", {
    list <- plotGene(bcb, genes = genes, return = "list")
    expect_is(list, "list")
    expect_true(
        lapply(list, function(x) is(x, "ggplot")) %>%
            unlist() %>%
            all()
    )
})
