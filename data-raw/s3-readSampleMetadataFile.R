library(tidyverse)
uploadDir <- system.file(
    file.path("extdata", "bcbio"),
    package = "bcbioRNASeq")
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
write_csv(meta, path = "~/Desktop/sample_metadata.csv")
