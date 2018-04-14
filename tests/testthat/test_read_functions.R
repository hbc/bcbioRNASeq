context("Read Functions")



# loadRNASeq ===================================================================
test_that("transformationLimit", {
    skip <- suppressWarnings(
        loadRNASeq(
            uploadDir = uploadDir,
            organism = "Mus musculus",
            transformationLimit = -Inf
        )
    )
    expect_identical(
        names(assays(skip)),
        c("counts", "tpm", "length", "normalized")
    )
})

test_that("User-defined sample metadata", {
    bcb <- suppressWarnings(loadRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        sampleMetadataFile = file.path(uploadDir, "sample_metadata.csv")
    ))
    expect_s4_class(bcb, "bcbioRNASeq")
    expect_identical(
        basename(metadata(bcb)[["sampleMetadataFile"]]),
        "sample_metadata.csv"
    )
})
