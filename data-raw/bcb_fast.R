## bcbioRNASeq fast mode example object.
## Updated 2022-03-07.
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(lobstr)
    library(basejump)
})
load_all()
## Restrict to 1 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(1e6, class = "object_size")
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
object <- bcbioRNASeq(uploadDir = uploadDir, fast = TRUE)
## Report the size of each slot in bytes.
lapply(coerceToList(object), obj_size)
stopifnot(
    isTRUE(obj_size(object) < limit),
    is(object, "bcbioRNASeq"),
    validObject(object)
)
## Upload this to S3 bucket.
assignAndSaveData(name = "bcb_fast", object = object, dir = "~")
