## bcbioRNASeq fast mode example object.
## Updated 2022-05-24.
## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(basejump)
})
## nolint end
load_all()
## Restrict to 1 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(1e6L, class = "object_size")
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
object <- bcbioRNASeq(uploadDir = uploadDir, fast = TRUE)
## Report the size of each slot in bytes.
stopifnot(
    isTRUE(object.size(object) < limit),
    is(object, "bcbioRNASeq"),
    validObject(object)
)
## Upload this to S3 bucket.
assignAndSaveData(name = "bcb_fast", object = object, dir = "~")
