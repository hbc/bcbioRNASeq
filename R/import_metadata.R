#' Import data from bcbio-nextgen run
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @rdname import
#' @description Import metadata
#'
#' @import dplyr
#' @import readr
#' @import stringr
#' @import tidyr
#' @importFrom utils write.csv
#'
#' @param bcbio \code{bcbio-nextgen} run object
#' @param grep Apply grep pattern matching to subset by description
#' @param save Save data frame
#' @param pool Pool lane split samples
#'
#' @return Metadata data frame
#' @export
import_metadata <- function(
    bcbio,
    grep = NULL,
    save = FALSE,
    pool = FALSE) {
    check_bcbio(bcbio)
    metadata <- list.files(bcbio$config_dir,
                           pattern = ".csv",
                           full.names = TRUE) %>%
        read_csv(col_types = cols()) %>%
        set_names_snake

    if (isTRUE(bcbio$lane_split) & !isTRUE(pool)) {
        # Lane splitting. This assumes the YAML descriptions won't match the
        # `_L00[1-4]` suffix. Therefore, it renames both the samplename and
        # description columns to match the bcbio server output. This workflow is
        # used by Harvard Biopolymers Facility. We can decide to either combine
        # counts at the server level using `cat` in bash, or we can run DESeq2
        # later by pooling the counts with `deseq_lane_pool()`. We may want to
        # deprecate this method in the future and simply combine counts at the
        # server level for all lane split runs.
        lane <- paste0("L", stringr::str_pad(1:4, 3, pad = "0"))

        metadata <- metadata %>%
            group_by_(.dots = "samplename") %>%
            expand_(~lane) %>%
            left_join(metadata, by = "samplename") %>%
            ungroup
        metadata$samplename <- paste(metadata$samplename, lane, sep = "_")

        # Check against the sample_dirs
        description_match <- identical(metadata$description,
                                       names(bcbio$sample_dirs))
        samplename_match <- identical(metadata$samplename,
                                      names(bcbio$sample_dirs))

        # Replace description for lanesplit samples
        if (!isTRUE(description_match) & isTRUE(samplename_match)) {
            metadata$description <- metadata$samplename
        }
    } else if (isTRUE(pool)) {
        metadata$description <- metadata$samplename
    }

    if (!is.null(grep)) {
        metadata <- metadata[grepl(grep, metadata$description), ]
    }

    metadata <- metadata %>% arrange_(.dots = "description")

    if (isTRUE(save)) {
        save(metadata, file = "data/metadata.rda")
        write.csv(metadata, file = "meta/metadata.csv")
    }

    return(metadata)
}
