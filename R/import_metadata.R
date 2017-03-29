#' Metadata for bcbio run
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#' @import stringr
#' @import tidyr
#' @importFrom basejump printTable
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

    # Lane splitting
    # Workflow used by Harvard Biopolymers Facility
    if (isTRUE(lane_split)) {
        lane <- paste0("L", stringr::str_pad(1:4, 3, pad = "0"))
        metadata <- metadata %>%
            dplyr::group_by_(.dots = "samplename") %>%
            tidyr::expand_(~lane) %>%
            dplyr::left_join(metadata, by = "samplename") %>%
            dplyr::ungroup() %>%
            dplyr::mutate_(.dots = set_names(
                list(quote(paste(samplename, lane, sep = "_"))),
                "samplename"
            ))
    }

    # Arrange the rows by description
    metadata <- metadata %>%
        dplyr::arrange_(.dots = "samplename") %>%
        set_rownames("samplename")

    if (isTRUE(save)) {
        save(metadata, file = "data/metadata.rda")
        write.csv(metadata, file = "results/metadata.csv")
        metadata[, c("samplename",
                     "description",
                     bcbio$intgroup)] %>%
            basejump::printTable(., caption = "Sample metadata")
    }

    return(metadata)
}
