context("tmm")

test_that("tmm", {
    for (object in list(
        bcbioRNASeq = object,
        matrix = assay(object)
    )) {
        expect_is(
            object = tmm(object),
            class = "matrix"
        )
    }
})
