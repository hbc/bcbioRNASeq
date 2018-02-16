#' Assert Checks
#'
#' @rdname assert
#' @name assert
#'
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than
#' @importFrom assertive assert_all_are_non_negative
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_are_intersecting_sets
#' @importFrom assertive assert_are_set_equal
#' @importFrom assertive assert_has_colnames
#' @importFrom assertive assert_has_dimnames
#' @importFrom assertive assert_has_rows
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_number
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_all_of
#' @importFrom assertive assert_is_an_integer
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_data.frame
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_matrix
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_numeric
#' @importFrom assertive assert_is_scalar
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl_df
#' @importFrom assertive assert_is_vector
#' @importFrom assertive is_a_string
#'
#' @importFrom basejump assert_formal_annotation_col
#' @importFrom basejump assert_formal_color_function
#' @importFrom basejump assert_formal_gene2symbol
#' @importFrom basejump assert_formal_header_level
#' @importFrom basejump assert_is_a_number_or_null
#' @importFrom basejump assert_is_a_string_or_null
#' @importFrom basejump assert_is_an_implicit_integer
#' @importFrom basejump assert_is_an_implicit_integer_or_null
#' @importFrom basejump assert_is_annotable
#' @importFrom basejump assert_is_character_or_null
#' @importFrom basejump assert_is_data.frame_or_null
#' @importFrom basejump assert_is_gene2symbol
#' @importFrom basejump assert_is_implicit_integer
#' @importFrom basejump assert_is_tx2gene
#'
#' @importFrom bcbioBase assert_formal_interesting_groups
#'
#' @inheritParams general
#'
#' @param severity How severe should the consequences of the assertion be?
#'   Either "`stop`", "`warning`", "`message`", or "`none`".
#'
#' @return Stop on error.
NULL



.assert_formal_scale_discrete <- function(  # nolint
    x,
    severity = "stop") {
    assert_is_any_of(
        x = x,
        classes = c("ScaleDiscrete", "NULL"),
        severity = severity)
}



.assert_formal_title <- function(  # nolint
    x,
    severity = "stop") {
    assert_is_any_of(
        x = x,
        classes = c("character", "logical", "NULL"),
        severity = severity)
}
