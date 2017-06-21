#' Import RNA-seq counts
#'
#' Import RNA-seq counts using [tximport()]. Currently supports
#' [salmon](https://combine-lab.github.io/salmon/) (**recommended**) and
#' [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/).
#'
#' @rdname tximport
#' @keywords internal
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param sample_dirs Sample directories to import.
#' @param tx2gene Transcript to gene annotations.
#'
#' @return Counts saved in [tximport] list object.
#'
#' @seealso
#' - [tximport::tximport()].
.tximport <- function(sample_dirs, tx2gene) {
    # Check for count output format, by using the first sample directory
    subdirs <- list.dirs(sample_dirs[1], full.names = FALSE, recursive = FALSE)
    if ("salmon" %in% subdirs) {
        type <- "salmon"
    } else if ("sailfish" %in% subdirs) {
        type <- "sailfish"
    } else {
        stop("Unsupported counts output format")
    }

    # Locate `quant.sf` file for salmon or sailfish output
    if (type %in% c("salmon", "sailfish")) {
        sample_files <- list.files(
            file.path(sample_dirs, type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE)
    }

    # Assign descriptions to sample files
    names(sample_files) <- names(sample_dirs)

    # Begin loading of selected counts
    message(paste("Reading", type, "counts"))

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    if (type %in% c("salmon", "sailfish")) {
        countsFromAbundance <- "lengthScaledTPM"
    } else {
        countsFromAbundance <- "no"
    }

    tximport(files = sample_files,
             type = type,
             tx2gene = as.data.frame(tx2gene),
             importer = read_tsv,
             countsFromAbundance = countsFromAbundance)
}



#' Read transcript to gene (tx2gene) annotation file
#'
#' @rdname tx2gene
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.tx2gene <- function(project_dir, genome_build = NULL) {
    file_path <- file.path(project_dir, "tx2gene.csv")
    if (file.exists(file_path)) {
        # Read bcbio-nextgene tx2gene file
        .read_file(file_path, col_names = c("enstxp", "ensgene")) %>%
            as.data.frame %>%
            # [TODO] Add S4 DataFrame method support for arrange in basejump
            arrange(!!sym("enstxp")) %>%
            DataFrame %>%
            set_rownames(.[["enstxp"]])
    } else {
        # Load prebuilt annotables tx2gene table
        if (is.null(genome_build)) {
            stop("Genome build required for annotable tx2gene")
        }
        .annotable(genome_build, format = "tx2gene")
    }
}
