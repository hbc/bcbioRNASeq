test_that("Facet wrapped", {
    p <- plotCounts(
        object = object,
        genes = geneNames,
        interestingGroups = "sampleName",
        style = "facet"
    )
    expect_s3_class(p, "ggplot")
})

test_that("Wide format", {
    p <- plotCounts(
        object = object,
        genes = geneNames,
        style = "wide"
    )
    expect_s3_class(p, "ggplot")
})
