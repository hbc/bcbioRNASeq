library(tidyverse)
uploadDir <- system.file(
    file.path("extdata", "bcbio"),
    package = "bcbioRNASeq")
bcb <- loadRNASeq(uploadDir)
meta <- sampleMetadata(bcb)
write_csv(meta, path = "~/Desktop/sample_metadata.csv")
