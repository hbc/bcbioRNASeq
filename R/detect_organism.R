#' Detect the organism from the genome build name
#'
#' @param genome_build Genome build
#'
#' @return Organism string
#' @export
#'
#' @examples
#' detect_organism("hg19")
#' detect_organism("mm10")
#' detect_organism("WBcel235")
#' detect_organism("BDGP6")
detect_organism <- function(genome_build) {
    if (str_detect(genome_build, "^(hg|GRCh)\\d+")) {
        organism <- "hsapiens"
    } else if (str_detect(genome_build, "^mm\\d+")) {
        organism <- "mmusculus"
    } else if (str_detect(genome_build, "^WBcel\\d+")) {
        organism <- "celegans"
    } else if (str_detect(genome_build, "^BDGP\\d+")) {
        organism <- "dmelanogaster"
    } else if (str_detect(genome_build, "^Zv\\d+")) {
        organism <- "drerio"
    } else if (str_detect(genome_build, "^ASM\\d+")) {
        organism <- "spombe"
    } else if (str_detect(genome_build, "^(MB|MG)\\d+")) {
        organism <- "ecoli"
    } else {
        message("Unknown genome build")
        organism <- NULL
    }
    return(organism)
}
