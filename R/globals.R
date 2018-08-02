globalVariables(".")

packageVersion <- packageVersion("bcbioRNASeq")
legacyMetricsCols <- c("name", "x53Bias")

# v0.2.6: tpm and length are now optional, since we're supporting featureCounts
requiredAssays <- "counts"

# v0.2.6: added STAR and HISAT2 support
tximportCallers <- c("salmon", "kallisto", "sailfish")
featureCountsCallers <- c("star", "hisat2")
validCallers <- c(tximportCallers, featureCountsCallers)
