library(tidyverse)
uploadDir <- system.file(
    file.path("extdata", "bcbio"),
    package = "bcbioRNASeq")
bcb <- loadRNASeq(uploadDir)
meta <- sampleMetadata(bcb) %>%
    # Add required `fileName` column
    mutate(fileName = paste0(sampleID, ".fastq.gz")) %>%
    select(fileName, everything())
write_csv(meta, path = "~/Desktop/sample_metadata.csv")
