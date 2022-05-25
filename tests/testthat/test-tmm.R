test_that("tmm", {
    for (object in list(
        "bcbioRNASeq" = object,
        "matrix" = assay(object)
    )) {
        expect_type(
            object = tmm(object),
            type = "double"
        )
    }
})
