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
#' @param intgroup Internal sample groups used for plot coloring.
#' @param lane_split Whether samples were split across flow cell lanes.
#'
#' @return Metadata data frame
#' @export
#'
#' @examples
#' \dontrun{
#' import_metadata(bcbio, intgroup = "treatment")
#' }
import_metadata <- function(
    bcbio,
    intgroup,
    lane_split = FALSE) {
    metadata <- list.files(bcbio$config_dir,
                           pattern = ".csv",
                           full.names = TRUE) %>%
        readr::read_csv(., col_types = readr::cols()) %>%
        basejump::setNamesSnake(.)

    # Set intgroup, if specified
    if (intgroup %in% colnames(metadata)) {
        metadata <- dplyr::mutate_(
            metadata,
            .dots = set_names(list(intgroup), "intgroup")
        )
    } else {
        stop("Required intgroup is not present in the config metadata.")
    }

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
                "description"
            ))
    }

    # Arrange the rows by description
    metadata <- metadata %>%
        dplyr::arrange_(.dots = "samplename") %>%
        set_rownames("samplename")

    return(metadata)
}
