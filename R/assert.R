#' Assert Checks
#'
#' @rdname assert
#' @name assert
#'
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than
#' @importFrom assertive assert_all_are_non_empty_character
#' @importFrom assertive assert_all_are_non_negative
#' @importFrom assertive assert_all_are_positive
#' @importFrom assertive assert_are_disjoint_sets
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
#' @importFrom assertive assert_is_factor
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_matrix
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_numeric
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl_df
#' @importFrom assertive assert_is_vector
#' @importFrom assertive is_a_string
#' @importFrom assertive is_positive
#'
#' @importFrom basejump assertFormalAnnotationCol
#' @importFrom basejump assertIsHexColorFunctionOrNULL
#' @importFrom basejump assertFormalGene2symbol
#' @importFrom basejump assertIsAHeaderLevel
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAnImplicitIntegerOrNULL
#' @importFrom basejump assertIsAnnotable
#' @importFrom basejump assertIsCharacterOrNULL
#' @importFrom basejump assertIsImplicitInteger
#' @importFrom basejump assertIsGene2symbol
#' @importFrom basejump assertIsColorScaleDiscreteOrNULL
#' @importFrom basejump assertIsFillScaleDiscreteOrNULL
#' @importFrom basejump assertIsTx2gene
#'
#' @importFrom bcbioBase assertFormalInterestingGroups
#'
#' @inheritParams general
#'
#' @param severity How severe should the consequences of the assertion be?
#'   Either "`stop`", "`warning`", "`message`", or "`none`".
#'
#' @return Stop on error.
NULL
