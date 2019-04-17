context("plotCountsPerGene")

geom <- methodFormals(
    f = "plotCountsPerGene",
    signature = "bcbioRNASeq"
) %>%
    .[["geom"]] %>%
    eval()

with_parameters_test_that(
    "plotCountsPerGene", {
        x <- plotCountsPerGene(object, geom = geom)
        expect_s3_class(x, "ggplot")
        x <- plotCountsPerGene(
            object = object,
            normalized = "vst",
            interestingGroups = "sampleName",
            title = NULL
        )
        expect_s3_class(x, "ggplot")
    },
    geom = geom
)
