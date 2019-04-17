context("tmm")

with_parameters_test_that(
    "tmm", {
        expect_is(
            object = tmm(object),
            class = "matrix"
        )
    },
    object = list(
        bcbioRNASeq = object,
        matrix = assay(object)
    )
)
