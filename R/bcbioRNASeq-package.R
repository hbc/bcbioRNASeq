#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#'
#' @import methods S4Vectors
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importFrom SummarizedExperiment assay assayNames assays colData rowData
#'   rowRanges
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than
#' @importFrom assertive assert_all_are_in_left_open_range
#' @importFrom assertive assert_all_are_matching_regex
#' @importFrom assertive assert_all_are_non_empty_character
#' @importFrom assertive assert_all_are_non_negative
#' @importFrom assertive assert_all_are_positive
#' @importFrom assertive assert_are_disjoint_sets
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_are_intersecting_sets
#' @importFrom assertive assert_are_set_equal
#' @importFrom assertive assert_has_colnames
#' @importFrom assertive assert_has_dimnames
#' @importFrom assertive assert_has_dims
#' @importFrom assertive assert_has_no_duplicates
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
#' @importFrom assertive assert_is_formula
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_matrix
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_numeric
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl_df
#' @importFrom assertive assert_is_vector
#' @importFrom assertive has_dims
#' @importFrom assertive has_names
#' @importFrom assertive is_a_string
#' @importFrom assertive is_positive
#' @importFrom basejump assertAreGeneAnnotations
#' @importFrom basejump assertFormalAnnotationCol
#' @importFrom basejump assertIsHexColorFunctionOrNULL
#' @importFrom basejump assertFormalGene2symbol
#' @importFrom basejump assertIsAHeaderLevel
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAnImplicitIntegerOrNULL
#' @importFrom basejump assertIsCharacterOrNULL
#' @importFrom basejump assertIsImplicitInteger
#' @importFrom basejump assertIsGene2symbol
#' @importFrom basejump assertIsColorScaleDiscreteOrNULL
#' @importFrom basejump assertIsFillScaleDiscreteOrNULL
#' @importFrom basejump assertIsTx2gene
#' @importFrom bcbioBase assertFormalInterestingGroups
#' @importFrom rlang .data abort inform warn
#' @importFrom utils globalVariables packageVersion
NULL



globalVariables(".")
packageVersion <- packageVersion("bcbioRNASeq")
lanePattern <- "_L(\\d{3})"
metadataPriorityCols <- c("sampleID", "description", "sampleName")
legacyMetricsCols <- c(metadataPriorityCols, "name", "x53Bias")
updateMsg <- "Run `updateObject()` to update your object"
validCallers <- c("salmon", "kallisto", "sailfish")
requiredAssays <- c("raw", "tpm", "length")
