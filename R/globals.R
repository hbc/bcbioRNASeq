globalVariables(".")
legacyMetricsCols <- c("name", "x53Bias")
packageVersion <- packageVersion("bcbioRNASeq")
# v0.2.6: tpm and length are now optional, since we're supporting featureCounts
requiredAssays <- "counts"
# v0.2.6: added STAR and HISAT2 support
validCallers <- c("salmon", "kallisto", "sailfish", "star", "hisat2")
