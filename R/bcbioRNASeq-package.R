#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#' @docType package
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_greater_than
#'   assert_all_are_in_left_open_range assert_all_are_in_range
#'   assert_all_are_non_negative assert_all_are_positive is_positive
#' @importFrom assertive.properties assert_has_colnames assert_has_dimnames
#'   assert_has_dims assert_has_no_duplicates assert_has_rows
#'   assert_is_non_empty assert_is_of_length assert_is_vector has_dims has_names
#' @importFrom assertive.sets assert_are_disjoint_sets
#'   assert_are_intersecting_sets assert_are_set_equal assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_empty_character
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_data.frame assert_is_factor assert_is_formula
#'   assert_is_list assert_is_matrix assert_is_numeric assert_is_tbl_df
#'   is_a_string
#' @importFrom basejump assertAreGeneAnnotations assertIsHexColorFunctionOrNULL
#'   assertFormalGene2symbol assertIsAHeaderLevel assertIsAStringOrNULL
#'   assertIsAnImplicitInteger assertIsAnImplicitIntegerOrNULL
#'   assertIsImplicitInteger assertIsGene2symbol
#'   assertIsColorScaleDiscreteOrNULL assertIsFillScaleDiscreteOrNULL
#'   assertIsTx2gene camel convertGenesToSymbols detectOrganism emptyRanges
#'   fixNA gene2symbol hasRownames initializeDirectory interestingGroups
#'   interestingGroups<- makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment markdownHeader markdownList markdownPlotlist
#'   plotHeatmap prepareTemplate readFileByExtension readYAML sampleData
#'   sanitizeRowData sanitizeSampleData snake stripTranscriptVersions
#'   uniteInterestingGroups matchInterestingGroups
#' @importFrom bcbioBase bcbio_geom_abline bcbio_geom_label
#'   bcbio_geom_label_repel copyToDropbox readDataVersions readLog
#'   readProgramVersions readSampleData readTx2gene readYAMLSampleData
#'   readYAMLSampleMetrics sampleDirs
#' @importFrom BiocGenerics cbind colSums density design
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom DEGreport degCovariates significants
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromMatrix
#'   DESeqDataSetFromTximport DESeqTransform estimateDispersions
#'   estimateSizeFactors results rlog varianceStabilizingTransformation
#' @importFrom dplyr arrange bind_cols desc everything filter group_by left_join
#'   mutate mutate_all mutate_if pull rename row_number select select_if
#'   starts_with ungroup
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom ggplot2 aes annotation_logticks coord_fixed coord_flip
#'   element_blank element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_density geom_hline geom_jitter geom_point geom_polygon geom_ribbon
#'   geom_smooth geom_vline ggplot ggtitle guides labs position_jitterdodge
#'   scale_color_manual scale_x_continuous scale_y_continuous stat_summary theme
#'   xlab ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom matrixStats colMedians
#' @importFrom methods .hasSlot as as<- is new setAs show slot slot<-
#'   validObject
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom reshape2 melt
#' @importFrom rlang := !! !!! sym syms
#' @importFrom S4Vectors as.data.frame complete.cases head mcols mcols<-
#'   metadata na.omit
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_match str_trunc
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#'   colData colData<- metadata<- rowData rowRanges SummarizedExperiment
#' @importFrom tibble as_tibble column_to_rownames glimpse remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom vsn meanSdPlot
NULL
