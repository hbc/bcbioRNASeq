globalVariables(".")
featureCountsAssays <- requiredAssays
featureCountsCallers <- c("star", "hisat2")
lanePattern <- basejump::lanePattern
legacyMetricsCols <- c("name", "x53Bias")
metadataBlacklist <- bcbioBase::metadataBlacklist
packageVersion <- packageVersion("bcbioRNASeq")
projectDirPattern <- bcbioBase::projectDirPattern
requiredAssays <- "counts"
separatorBar <- basejump::separatorBar
tximportAssays <- c("counts", "length", "tpm")
tximportCallers <- c("salmon", "kallisto", "sailfish")
updateMessage <- basejump::updateMessage
validCallers <- c(tximportCallers, featureCountsCallers)
validLevels <- c("genes", "transcripts")
