context("aggregateReplicates")

test_that("matrix", {
    genes <- c("ENSMUSG00000000001",
               "ENSMUSG00000000003",
               "ENSMUSG00000000028",
               "ENSMUSG00000000031",
               "ENSMUSG00000000037")
    # This function works by detecting lane split replicates (e.g. `_L001`)
    counts <- data.frame(
        sample1_L001 = c(0, 1, 2, 3, 4),
        sample1_L002 = c(0, 1, 2, 3, 4),
        sample1_L003 = c(0, 1, 2, 3, 4),
        sample1_L004 = c(0, 1, 2, 3, 4),
        sample2_L001 = c(4, 3, 2, 1, 0),
        sample2_L002 = c(4, 3, 2, 1, 0),
        sample2_L003 = c(4, 3, 2, 1, 0),
        sample2_L004 = c(4, 3, 2, 1, 0),
        row.names = genes
    ) %>%
        as.matrix()
    aggregated <- data.frame(
        sample1 = c(0, 4, 8, 12, 16),
        sample2 = c(16, 12, 8, 4, 0),
        row.names = genes
    ) %>%
        as.matrix()
    expect_equal(
        aggregateReplicates(counts),
        aggregated
    )
})
