## bcbioRNASeq example object.
## Updated 2022-05-24.
## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(basejump)
})
## nolint end
load_all()
## Restrict to 2 MB.
limit <- structure(2e6L, class = "object_size")
uploadDir <- system.file("extdata", "bcbio", package = "bcbioRNASeq")
object <- bcbioRNASeq(
    uploadDir = uploadDir,
    level = "genes",
    caller = "salmon",
    interestingGroups = c("treatment", "day"),
    organism = "Mus musculus",
    ensemblRelease = 90L
)
## DESeq2 doesn't like spaces in design formula factors.
object[["treatment"]] <- snakeCase(object[["treatment"]])
## Report the size of each slot in bytes.
stopifnot(
    isTRUE(object.size(object) < limit),
    is(object, "bcbioRNASeq"),
    validObject(object)
)
bcb <- object
use_data(bcb, overwrite = TRUE, compress = "xz")
