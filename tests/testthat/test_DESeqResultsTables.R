lfcThreshold <- 0.25

# export
# markdown

# FIXME This needs an update.
# test_that("resultsTables : DESeqAnalysis", {
#     output <- capture.output(
#         resultsTables(
#             object = res_small,
#             counts = dds_small,
#             lfcThreshold = lfcThreshold,
#             summary = TRUE,
#             headerLevel = 2L,
#             write = TRUE,
#             dir = "resultsTables"
#         ),
#         type = "output"
#     ) %>%
#         # Just evaluate the R Markdown output
#         head(28L)
#     expect_identical(
#         output,
#         c(
#             "",
#             "",
#             "## treatment folic acid vs control",
#             "",
#             "",
#             "",
#             "### Summary statistics",
#             "",
#             "",
#             "- 500 genes in counts matrix",
#             "- Base mean > 0: 500 genes (non-zero)",
#             "- Base mean > 1: 500 genes",
#             "- Alpha: 0.1",
#             "- LFC threshold: 0.25",
#             "- DEG pass alpha: 254 genes",
#             "- DEG LFC up: 115 genes",
#             "- DEG LFC down: 139 genes",
#             "",
#             "",
#             "",
#             "### Results tables",
#             "",
#             "",
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_all.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_all.csv.gz"
#                 ),
#                 "): ",
#                 "All genes, sorted by Ensembl identifier."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg.csv.gz"
#                 ),
#                 "): ",
#                 "Genes that pass the alpha (FDR) cutoff."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg_lfc_up.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg_lfc_up.csv.gz"
#                 ),
#                 "): ",
#                 "Upregulated DEG; positive log2 fold change."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg_lfc_down.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg_lfc_down.csv.gz"
#                 ),
#                 "): ",
#                 "Downregulated DEG; negative log2 fold change."
#             ),
#             ""
#         )
#     )
# })
#
# if (file.exists("token.rds")) {
#     test_that("resultsTables : DESeqAnalysis : Dropbox mode", {
#         resTbl <- resultsTables(
#             results = res_small,
#             counts = dds_small,
#             lfcThreshold = lfcThreshold,
#             summary = FALSE,
#             write = TRUE,
#             dir = "resultsTables",
#             dropboxDir = file.path("bcbioRNASeq_examples", "resultsTables"),
#             rdsToken = "token.rds"
#         )
#         expect_true("dropboxFiles" %in% names(resTbl))
#         # Check for Dropbox URLs
#         expect_true(all(vapply(
#             X = resTbl[["dropboxFiles"]],
#             FUN = function(file) {
#                 grepl("^https://www.dropbox.com/s/", file[["url"]])
#             },
#             FUN.VALUE = logical(1L)
#         )))
#         # Now check the Markdown code
#         output <- capture.output(.markdownResultsTables(resTbl))
#         # The function currently returns links to 4 DEG files
#         expect_identical(
#             length(which(grepl("https://www.dropbox.com/s/", output))),
#             4L
#         )
#     })
# }

