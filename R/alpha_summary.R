#' Print summary statistics of alpha level cutoffs
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param dds \code{DESeqDataSet}
#' @param alpha alpha Numeric vector of alpha levels to check
#'
#' @return Printed summary
#' @export
alpha_summary <- function(
    dds,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6)) {
    name <- deparse(substitute(dds))
    message(name)
    lapply(seq_along(alpha), function(a) {
        writeLines(paste0(name, ": alpha = ", alpha[a]))
        results(dds, alpha = alpha[a]) %>% summary
    }) %>% invisible
}
