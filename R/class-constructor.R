.bcbioRnaDataSet <- function(
    analysis,
    counts,
    alt_counts = NULL,
    groups_of_interest,
    custom_metadata_file = NULL) {
    # [TODO] Generate metadata from file
    custom_metadata <- NULL
    new("bcbioRnaDataSet",
        analysis = analysis,
        date = Sys.Date(),  # auto
        counts = counts,
        alt_counts = alt_counts,
        groups_of_interest = groups_of_interest,
        custom_metadata_file = custom_metadata_file,
        custom_metadata = custom_metadata,
        wd = getwd(),  # auto
        hpc = detect_hpc(),  # auto
        session_info = sessionInfo()  # auto
    )
}
