.new.bcbioRNASeq <- function(  # nolint
    assays,
    rowRanges,
    colData,
    metadata,
    transgeneNames,
    spikeNames
) {
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
    assert_is_all_of(rse, "RangedSummarizedExperiment")
    new("bcbioRNASeq", rse)
}
