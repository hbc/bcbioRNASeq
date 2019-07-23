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

normalizedCounts <- c("vst", "rlog", "tmm", "rle", "tpm", "fpkm")

legacyMetricsCols <- c("name", "x53Bias")

trans <- c("log2", "log10")

Rle <- structure("Rle", package = "S4Vectors")  # nolint
