context("plotCountsPerFeature")

geom <- methodFormals(
    f = "plotCountsPerFeature",
    signature = "bcbioRNASeq"
) %>%
    .[["geom"]] %>%
    eval()

with_parameters_test_that(
    "bcbioRNASeq", {
        x <- plotCountsPerFeature(object, geom = geom)
        expect_s3_class(x, "ggplot")
        x <- plotCountsPerFeature(
            object = object,
            normalized = "vst",
            interestingGroups = "sampleName",
            title = NULL
        )
        expect_s3_class(x, "ggplot")
    },
    geom = geom
)



context("plotCountsPerGene")

test_that("bcbioRNASeq", {
    x <- plotCountsPerGene(object)
    expect_s3_class(x, "ggplot")
})
