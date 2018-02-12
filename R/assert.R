#' Assert Checks
#'
#' @rdname assert
#' @name assert
#'
#' @inheritParams general
NULL



#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_data.frame
#' @importFrom assertive assert_is_matrix
#' @importFrom assertive assert_is_subset



#' @importFrom bcbioBase checkGene2symbol
.assert_gene2symbol <- function(object, genes, gene2symbol) {
    assert_is_any_of(genes, c("character", "NULL"))
    if (is.character(genes)) {
        assert_is_subset(genes, rownames(object))
    }
    assert_is_any_of(gene2symbol, c("data.frame", "NULL"))
    if (is.data.frame(gene2symbol)) {
        # TODO Rename to `assert_is_gene2symbol`
        checkGene2symbol(gene2symbol)
        assert_is_subset(rownames(object), rownames(gene2symbol))
    }
}
