# Initialize interesting groups
.prepareInterestingGroups <- function(
    object,
    interestingGroups
) {
    if (missing(interestingGroups)) {
        interestingGroups <- basejump::interestingGroups(object)
        if (is.null(interestingGroups)) {
            interestingGroups <- "sampleName"  # nocov
        }
    } else if (is.null(interestingGroups)) {
        interestingGroups <- "sampleName"
    } else {
        interestingGroups(object) <- interestingGroups
    }
    assert_is_character(interestingGroups)
    assert_is_non_empty(interestingGroups)
    interestingGroups
}
