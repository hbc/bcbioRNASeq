#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#' @keywords internal
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importFrom BiocGenerics colSums density design
#' @importFrom DEGreport degCovariates
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromTximport DESeqTransform
#'   estimateSizeFactors results rlog varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors as.data.frame head mcols metadata na.omit
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#'   colData colData<- rowData rowRanges SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols detectOrganism emptyRanges
#'   fixNA initializeDirectory makeGRangesFromEnsembl makeGRangesFromGFF
#'   makeNames makeTx2geneFromGFF markdownHeader markdownList markdownPlotlist
#'   readYAML sanitizeRowData sanitizeSampleData snake
#' @importFrom bcbioBase copyToDropbox flatFiles gene2symbol interestingGroups
#'   interestingGroups<- plotHeatmap prepareSummarizedExperiment prepareTemplate
#'   readDataVersions readLog readProgramVersions readSampleData readTx2gene
#'   sampleData sampleDirs sampleYAMLMetadata sampleYAMLMetrics
#'   uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols desc filter group_by left_join mutate
#'   mutate_all mutate_if pull rename select_if
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits facet_wrap geom_bar
#'   geom_boxplot geom_density geom_hline geom_jitter geom_point geom_polygon
#'   geom_ribbon geom_smooth ggplot ggtitle guides labs position_jitterdodge
#'   scale_color_manual scale_fill_manual scale_x_continuous scale_x_log10
#'   stat_summary theme xlab ylim
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom matrixStats colMedians
#' @importFrom methods .hasSlot as as<- is new show slot slot<- validObject
#' @importFrom parallel mclapply mcmapply
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom reshape2 melt
#' @importFrom rlang !! !!! sym syms
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble as_tibble column_to_rownames glimpse remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom viridis inferno scale_color_viridis scale_fill_viridis viridis
#' @importFrom vsn meanSdPlot
#'
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
#'
#' @importFrom basejump assertAreGeneAnnotations
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
#' @importFrom basejump hasRownames
#'
#' @importFrom bcbioBase assertFormalAnnotationCol
#' @importFrom bcbioBase assertFormalInterestingGroups
NULL
