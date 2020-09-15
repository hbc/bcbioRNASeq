## bcbioRNASeq fast mode example object.
## Updated 2020-09-14.

library(usethis)
library(pryr)

## Restrict to 1 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(1e6, class = "object_size")

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")

bcb <- bcbioRNASeq(
    uploadDir = uploadDir,
    fast = TRUE
)

object_size(bcb)
## 166 kB

## Report the size of each slot in bytes.
lapply(coerceS4ToList(bcb), object_size)
object_size(bcb)
stopifnot(object_size(bcb) < limit)

## Check that object is valid.
stopifnot(is(bcb, "bcbioRNASeq"))
validObject(bcb)

assignAndSaveData(name = "bcb_fast", object = bcb, dir = "~")
