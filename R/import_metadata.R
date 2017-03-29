#' Metadata for bcbio run
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#' @import stringr
#' @import tidyr
#'
#' @param bcbio bcbio run object
#' @param lane_split Whether samples were split across flow cell lanes
#' @param save Save data frame
#'
#' @return Metadata data frame
#' @export
#'
#' @examples
#' \dontrun{
#' import_metadata(bcbio)
#' }
import_metadata <- function(
    bcbio,
    lane_split = NULL,
    save = FALSE) {
    # Automatically detect if lanes are split, from bcbio object
    if (is.null(lane_split)) {
        lane_split <- bcbio$lane_split
    }

    metadata <- list.files(bcbio$config_dir,
                           pattern = ".csv",
                           full.names = TRUE) %>%
        readr::read_csv(., col_types = readr::cols()) %>%
        basejump::setNamesSnake(.)

    # Lane splitting This assumes the YAML descriptions won't match the
    # `_L00[1-4]` suffix. Therefore, it renames both the samplename and
    # description columns to match the bcbio server output. This workflow is
    # used by Harvard Biopolymers Facility. We can decide to either combine
    # counts at the server level using `cat` in bash, or we can run DESeq2 later
    # by pooling the counts with `deseq_lane_pool()`. We may want to deprecate
    # this method in the future and simply combine counts at the server level
    # for all lane split runs.
    if (isTRUE(lane_split)) {
        lane <- paste0("L", stringr::str_pad(1:4, 3, pad = "0"))
        metadata <- metadata %>%
            dplyr::group_by_(.dots = "samplename") %>%
            tidyr::expand_(~lane) %>%
            dplyr::left_join(metadata, by = "samplename") %>%
            dplyr::ungroup()
        metadata$samplename <- paste(metadata$samplename, lane, sep = "_")
        metadata$description <- metadata$samplename
    }

    # Ensure that rownames are set
    rownames(metadata) <- metadata$description

    if (isTRUE(save)) {
        save(metadata, file = "data/metadata.rda")
        write.csv(metadata, file = "results/metadata.csv")
    }

    return(metadata)
}
