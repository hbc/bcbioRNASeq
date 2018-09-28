globalVariables(".")

packageVersion <- packageVersion("bcbioRNASeq")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioRNASeqCacheURL
bcbioRNASeqCacheURL <- paste0(
    "http://bcbiobase.seq.cloud/",
    "v", packageVersion$major, ".", packageVersion$minor
)

lanePattern <- basejump::lanePattern
separatorBar <- basejump::separator()
updateMessage <- basejump::updateMessage

metadataBlacklist <- bcbioBase::metadataBlacklist
projectDirPattern <- bcbioBase::projectDirPattern

validLevels <- c("genes", "transcripts")

requiredAssays <- "counts"
tximportAssays <- c("counts", "length", "tpm")
featureCountsAssays <- requiredAssays

tximportCallers <- c("salmon", "kallisto", "sailfish")
featureCountsCallers <- c("star", "hisat2")
validCallers <- c(tximportCallers, featureCountsCallers)

legacyMetricsCols <- c("name", "x53Bias")
