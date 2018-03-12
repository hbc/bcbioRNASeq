library(tidyverse)
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
bcb <- loadRNASeq(uploadDir)
meta <- sampleMetadata(bcb) %>%
    mutate(
        # Add required `fileName` column
        fileName = paste0(sampleID, ".fastq.gz"),
        # Ensure `sampleID` and `sampleName` aren't defined
        sampleID = NULL,
        sampleName = NULL
    ) %>%
    select(fileName, everything())
write_csv(meta, path = "tests/testthat/sample_metadata.csv")
