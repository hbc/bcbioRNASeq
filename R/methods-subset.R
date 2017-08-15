#rdname bcbioRNADataSet
setMethod("[", c("bcbioRNADataSet", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE)
          {
              if (missing(i)) {
                  i <- 1L:nrow(x)
              }
              if (missing(j)) {
                  i <- 1L:ncol(x)
              }
              dots <- list(...)
              if (is.null(dots[["maxSamples"]]))
                  maxSamples <- 50L
              else maxSamples <- dots[["maxSamples"]]

              # Subset SE object ====
              se <- SummarizedExperiment(assays =
                                       SimpleList(counts(x)),
                                       colData = colData(x),
                                       metadata = metadata(x))

              tmp <- se[i, j, drop = drop]
              tmpGenes <- row.names(tmp)
              tmpData <- colData(tmp) %>% as.data.frame
              tmpTxi <- bcbio(x, "tximport")

              # Subset tximport data ====
              txi <- S4Vectors::SimpleList(
                  abundance = tmpTxi[["abundance"]] %>%
                      .[tmpGenes, tmpData[["sampleID"]]],
                  counts = tmpTxi[["counts"]] %>%
                      .[tmpGenes, tmpData[["sampleID"]]],
                  length = tmpTxi[["length"]] %>%
                      .[tmpGenes, tmpData[["sampleID"]]],
                  countsFromAbundance = tmpTxi[["countsFromAbundance"]]
              )
              rawCounts <- txi[["counts"]]
              tmm <- .tmm(rawCounts)
              tpm <- txi[["abundance"]]

              dds <- DESeqDataSetFromTximport(
                  txi = txi,
                  colData = tmpData,
                  design = formula(~1L)) %>%
                  DESeq
              normalizedCounts <- counts(dds, normalized = TRUE)

              # rlog & variance ====
              if (nrow(tmpData) > maxSamples){
                  message("Data to big, skipping vst/rlog")
                  rlog <- DESeqTransform(
                      SummarizedExperiment(assays = log2(tmm + 1L),
                                           colData = colData(dds)))
                  vst <- DESeqTransform(
                      SummarizedExperiment(assays = log2(tmm + 1L),
                                           colData = colData(dds)))
              }else{
                  message("Performing rlog transformation")
                  rlog <- rlog(dds)
                  message("Performing variance stabilizing transformation")
                  vst <- varianceStabilizingTransformation(dds)
              }

              if (is.matrix(bcbio(x, "featureCounts"))) {
                  tmpFC <- bcbio(x, "featureCounts") %>%
                      .[tmpGenes, tmpData[["sampleID"]]]
              }else{
                  tmpFC <- bcbio(x, "featureCounts")
              }

              # Subset Metrics ====
              tmpMetrics <- metadata(x)[["metrics"]] %>%
                  .[.[["sampleID"]] %in% tmpData[["sampleID"]], ]
              metadata(tmp)[["metrics"]] <- tmpMetrics

              assays(tmp) <- SimpleList(
                          raw = rawCounts,
                          normalized = normalizedCounts,
                          tpm = tpm,
                          tmm = tmm,
                          rlog = rlog,
                          vst = vst)

              # bcbioRNADataSet ====
              bcb <- new("bcbioRNADataSet", tmp)
              # Slot additional callers
              bcbio(bcb, "tximport") <- txi
              bcbio(bcb, "DESeqDataSet") <- dds
              if (is.matrix(tmpFC)) {
                  bcbio(bcb, "featureCounts") <- tmpFC
              }
              bcb

          })
