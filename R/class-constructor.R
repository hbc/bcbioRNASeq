.bcbioRnaDataSet <- function(
    analysis,
    counts,
    alt_counts = NULL,
    groups_of_interest,
    custom_metadata_file = NULL) {
    # [TODO] Generate metadata from file
    custom_metadata <- NULL
    bcb <- new("bcbioRnaDataSet",
        analysis = analysis,
        date = Sys.Date(),  # auto
        counts = counts,
        groups_of_interest = groups_of_interest,
        wd = getwd(),  # auto
        hpc = detect_hpc(),  # auto
        session_info = sessionInfo()  # auto
    )

    if (!is.null(alt_counts)) {
        bcb@alt_counts <- alt_counts
    } else if (!is.null(custom_metadata_file)) {
        bcb@custom_metadata_file <- custom_metadata_file
        # [TODO] Add workflow to process custom metadata
        custom_metadata <- DataFrame()
        bcb@custom_metadata <- custom_metadata
    }

    bcb
}
