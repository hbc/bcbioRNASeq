#' Detect the organism from the genome build name
#' @rdname internal-detect_organism
#' @param genome_build Genome build.
#' @return Organism string.
.detect_organism <- function(genome_build) {
    if (str_detect(genome_build, "^(hg|GRCh)\\d+")) {
        "hsapiens"
    } else if (str_detect(genome_build, "^mm\\d+")) {
        "mmusculus"
    } else if (str_detect(genome_build, "^WBcel\\d+")) {
        "celegans"
    } else if (str_detect(genome_build, "^BDGP\\d+")) {
        "dmelanogaster"
    } else if (str_detect(genome_build, "^Zv\\d+")) {
        "drerio"
    } else if (str_detect(genome_build, "^ASM\\d+")) {
        "spombe"
    } else if (str_detect(genome_build, "^(MB|MG)\\d+")) {
        "ecoli"
    } else {
        stop("Failed to detect organism")
    }
}
