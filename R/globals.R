globalVariables(".")

packageVersion <- packageVersion("bcbioRNASeq")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioRNASeqCacheURL
bcbioRNASeqCacheURL <- paste0(
    "http://bcbiornaseq.seq.cloud/",
    "v", packageVersion$major, ".", packageVersion$minor  # nolint
)

# Using internal list defined in basejump.
formalsList <- basejump:::formalsList

validLevels <- c("genes", "transcripts")

requiredAssays <- "counts"
tximportAssays <- c(requiredAssays, "tpm")
featureCountsAssays <- requiredAssays

tximportCallers <- c("salmon", "kallisto", "sailfish")
featureCountsCallers <- c("star", "hisat2")
validCallers <- c(tximportCallers, featureCountsCallers)

normalizedCounts <- c("vst", "rlog", "tmm", "rle", "tpm", "fpkm")

legacyMetricsCols <- c("name", "x53Bias")

Rle <- structure("Rle", package = "S4Vectors")  # nolint
