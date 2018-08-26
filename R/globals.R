globalVariables(".")

requiredAssays <- "counts"
featureCountsCallers <- c("star", "hisat2")
lanePattern <- basejump::lanePattern
legacyMetricsCols <- c("name", "x53Bias")
metadataBlacklist <- bcbioBase::metadataBlacklist
packageVersion <- packageVersion("bcbioRNASeq")
projectDirPattern <- bcbioBase::projectDirPattern
separatorBar <- basejump::separatorBar
tximportAssays <- c("counts", "length", "tpm")
tximportCallers <- c("salmon", "kallisto", "sailfish")
updateMessage <- basejump::updateMessage
validLevels <- c("genes", "transcripts")

featureCountsAssays <- requiredAssays
validCallers <- c(tximportCallers, featureCountsCallers)
