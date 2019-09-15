## FIXME Prefix non-exported globals with a ".".

globalVariables(".")

.version <- packageVersion("bcbioRNASeq")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioRNASeqTestsURL
bcbioRNASeqTestsURL <- paste0(
    "http://tests.acidgenomics.com/bcbioRNASeq/",
    "v", .version$major, ".", .version$minor  # nolint
)

validLevels <- c("genes", "transcripts")

requiredAssays <- "counts"
tximportAssays <- c(requiredAssays, "tpm")
featureCountsAssays <- requiredAssays

tximportCallers <- c("salmon", "kallisto", "sailfish")
featureCountsCallers <- c("star", "hisat2")
validCallers <- c(tximportCallers, featureCountsCallers)

## DESeqTransform.
dt <- c("vst", "rlog")

## "sf" denotes size-factor adjusted counts. This corresponds to
## `normalized = TRUE`, and is recommended by default for plots.
normalizedCounts <- c("tpm", "sf", "fpkm", dt, "tmm", "rle")

legacyMetricsCols <- c("name", "x53Bias")

trans <- c("log2", "log10")

Rle <- structure("Rle", package = "S4Vectors")  # nolint
