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
#' @importFrom BiocGenerics colSums density design
#' @importFrom DEGreport degCovariates significants
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromTximport DESeqTransform
#'   estimateSizeFactors results rlog varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors as.data.frame complete.cases head mcols mcols<-
#'   metadata na.omit
#' @importFrom SummarizedExperiment assay assay<- assayNames assays assays<-
#'   colData colData<- metadata<- rowData rowRanges SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols detectOrganism emptyRanges
#'   fixNA initializeDirectory makeGRangesFromEnsembl makeGRangesFromGFF
#'   makeNames makeTx2geneFromGFF markdownHeader markdownList markdownPlotlist
#'   readYAML sanitizeRowData sanitizeSampleData snake
#' @importFrom bcbioBase copyToDropbox flatFiles gene2symbol interestingGroups
#'   interestingGroups<- plotHeatmap prepareSummarizedExperiment prepareTemplate
#'   readDataVersions readLog readProgramVersions readSampleData readTx2gene
#'   readYAMLSampleData readYAMLSampleMetrics sampleData sampleDirs
#'   uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols desc filter group_by left_join mutate
#'   mutate_all mutate_if pull rename row_number select_if
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits facet_wrap geom_bar
#'   geom_boxplot geom_density geom_hline geom_jitter geom_point geom_polygon
#'   geom_ribbon geom_smooth geom_vline ggplot ggtitle guides labs
#'   position_jitterdodge scale_color_hue scale_color_manual scale_fill_hue
#'   scale_fill_manual scale_x_continuous scale_y_continuous stat_summary theme
#'   xlab ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom matrixStats colMedians
#' @importFrom methods .hasSlot as as<- is new show slot slot<- validObject
#' @importFrom parallel mclapply mcmapply
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom reshape2 melt
#' @importFrom rlang := !! !!! sym syms
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble as_tibble column_to_rownames glimpse remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom vsn meanSdPlot
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs
#' @importFrom assertive.files assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_greater_than
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#' @importFrom assertive.numbers assert_all_are_in_range
#' @importFrom assertive.numbers assert_all_are_non_negative
#' @importFrom assertive.numbers assert_all_are_positive
#' @importFrom assertive.numbers is_positive
#' @importFrom assertive.properties assert_has_colnames
#' @importFrom assertive.properties assert_has_dimnames
#' @importFrom assertive.properties assert_has_dims
#' @importFrom assertive.properties assert_has_no_duplicates
#' @importFrom assertive.properties assert_has_rows
#' @importFrom assertive.properties assert_is_non_empty
#' @importFrom assertive.properties assert_is_of_length
#' @importFrom assertive.properties assert_is_vector
#' @importFrom assertive.properties has_dims
#' @importFrom assertive.properties has_names
#' @importFrom assertive.sets assert_are_disjoint_sets
#' @importFrom assertive.sets assert_are_intersecting_sets
#' @importFrom assertive.sets assert_are_set_equal
#' @importFrom assertive.sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#' @importFrom assertive.strings assert_all_are_non_empty_character
#' @importFrom assertive.types assert_is_a_bool
#' @importFrom assertive.types assert_is_a_number
#' @importFrom assertive.types assert_is_a_string
#' @importFrom assertive.types assert_is_all_of
#' @importFrom assertive.types assert_is_an_integer
#' @importFrom assertive.types assert_is_any_of
#' @importFrom assertive.types assert_is_character
#' @importFrom assertive.types assert_is_data.frame
#' @importFrom assertive.types assert_is_factor
#' @importFrom assertive.types assert_is_formula
#' @importFrom assertive.types assert_is_list
#' @importFrom assertive.types assert_is_matrix
#' @importFrom assertive.types assert_is_numeric
#' @importFrom assertive.types assert_is_tbl_df
#' @importFrom assertive.types is_a_string
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
