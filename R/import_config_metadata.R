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
#' @param group Sample grouping. Used later in pipline for plot colors.
#' @param lane_split Whether samples were split across flow cell lanes.
#'
#' @return Metadata data frame
#' @export
#'
#' @examples
#' \dontrun{
#' import_config_metadata(bcbio, group = "treatment")
#' }
import_config_metadata <-
    function(project,
             group = "treatment",
             laneSplit = FALSE) {
        metadata <- list.files(project$configDir,
                               pattern = ".csv",
                               full.names = TRUE) %>%
            readr::read_csv(.,
                            col_types = readr::cols()) %>%
            basejump::setNamesSnake(.) %>%
            #` dplyr::rename_(.dots = c("name" = "samplename")) %>%
            dplyr::mutate_(.dots = set_names(list(group), "group"))

        # Lane splitting
        # Workflow used by Harvard Biopolymers Facility
        if (isTRUE(laneSplit)) {
            lane <- paste0("L", stringr::str_pad(1:4, 3, pad = "0"))
            metadata <- metadata %>%
                dplyr::group_by_(.dots = "samplename") %>%
                tidyr::expand_(~lane) %>%
                dplyr::left_join(metadata, by = "samplename") %>%
                dplyr::ungroup() %>%
                # NSE version
                # dplyr::mutate(description = paste(samplename, lane, sep = "_")) %>%
                dplyr::mutate_(.dots = set_names(
                    list(quote(paste(samplename, lane, sep = "_"))),
                    "description"
                ))
        }

        # Arrange the rows by description
        metadata <- metadata %>%
            dplyr::arrange_(.dots = "samplename") %>%
            basejump::setRownames(., "description")

        return(metadata)
    }
