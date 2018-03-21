#' @importFrom bcbioBase prepareSummarizedExperiment
.new.bcbioRNASeq <- function(  # nolint
    assays,
    rowRanges,
    colData,
    metadata,
    isSpike
) {
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
    assert_is_all_of(rse, "RangedSummarizedExperiment")
    new("bcbioRNASeq", rse)
}
