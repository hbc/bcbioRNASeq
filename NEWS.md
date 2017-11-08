# bcbioRNASeq 0.1.2 (2017-11-08)

- Updated package imports to match Bioconductor 3.6.
- Added support for interesting groups assignment with `interestingGroups<-`.
- Renamed `plotGeneHeatmap()` to simply `plotHeatmap()`.
- Added gender marker support for *Homo sapiens*.
- Improved support for multiple interesting groups in quality control plots. Now interestingGroups is defined as a column in the metrics data.frame that is used to specify the plot color/fill. This matches the convention in the bcbioSingleCell 0.0.22 update.
- Sample metadata columns are now consistently set as factors.



# bcbioRNASeq 0.1.1 (2017-10-26)

- Added support for coloring of multiple interesting groups in quality control plots.



# bcbioRNASeq 0.1.0 (2017-10-23)

- Updated version and author information to match the F1000 Research workflow.
- Added an `f1000v1` branch containing the reproducible code used to generate the figures in our workflow.
- Modified `plotMA()` to support vertical or horizontal layout return. Also added an argument to remove the color legend, which is typically not that informative.
- Added custom color palette support to the quality control functions.
- Upgrading from `bcbioRNADataSet` (< 0.1.0) to `bcbioRNASeq` class object is now possible using `as()` coercion method.
- Object oriented methods are now restricted to use `bcbioRNASeq` object. Legacy `bcbioRNADataSet` objects must be upgraded to `bcbioRNASeq` class.



# bcbioRNASeq 0.0.28 (2017-10-17)

- Added support for output of unstructured data inside `bcbioRNASeq` S4 object using `flatFiles()` function.
- Added `bcbioRNASeq` method support for `annotable()` generic.



# bcbioRNASeq 0.0.27 (2017-10-10)

- Renamed `bcbioRNADataSet` S4 class to `bcbioRNASeq`. This matches the naming conventions in the [bcbioSingleCell][] package.
- Renamed `loadRNASeqRun()` to simply `loadRNASeq()`.
- Switched `loadRNASeq()` from using S4 dispatch to a standard function.
- Added a parameter argument to `loadRNASeq()` that enables request of a specific Ensembl release version for gene annotations.
- Renamed `interestingGroup` argument in quality control functions to `interestingGroups` for better consistency.
- Improved handling of sample metrics in `plotPCACovariates()`.
- Added functional analysis R Markdown template.
- Offloaded some core functionality shared between [bcbioRNASeq][] and [bcbioSingleCell][] to the [basejump][] package. This included some code to handle sample metadata YAML and file loading. This helps provide a consistent experience across both packages.



# bcbioRNASeq 0.0.26 (2017-09-09)

- Renamed package from bcbioRnaseq to bcbioRNASeq.
- Improved website appearance.
- Added [viridis][] color palette support to quality control functions.
- Improved subset operations on `bcbioRNADataSet` object.
- Fixed setup chunk loading of `bcbioRNADataSet` in differential expression [RMarkdown][] template.



# bcbioRNASeq 0.0.25 (2017-08-11)

- Added S4 methods support for plots, allowing the user to use either
  `bcbioRNADataSet` or a metrics `data.frame` and manual `interesting_group`
  declaration for visualization.
- Migrated function and variable names from `snake_case` to `camelCase`.
- Offloaded small RNA functionality to a separate package named [bcbioSmallRNA][].



# bcbioRNASeq 0.0.24 (2017-07-13)

- Reworked [RMarkdown][] templates to improve YAML defaults and add more comments.
- Modified default path variables in `setup.R` to use `*_dir` instead of `*_out`.
- Updated NEWS file to use [Markdown][] syntax.



# bcbioRNASeq 0.0.23 (2017-07-03)

- Slotted `DESeqDataSet` using `design = formula(~1)` for quality control. This enables automatic generation of `rlog` and `vst` transformed counts.
- Documentation fixes and website updates.
- Renamed S4 class from `bcbioRnaDataSet` to `bcbioRNADataSet` (case sensitive).
- Adjusted the number of exported functions.



# bcbioRNASeq 0.0.22 (2017-06-21)

- Added [testthat][] checking with [lintr][].
- Initial setup of code coverage using [covr][].



# bcbioRNASeq 0.0.21 (2017-06-16)

- Prepared draft of [F1000][] workflow document.



# bcbioRNASeq 0.0.20 (2017-06-09)

- Added [Travis-CI][] support for automatic rendering of quality control report.



# bcbioRNASeq 0.0.19 (2017-06-07)

- `bcbioRnaDataSet` S4 definition updates.
- Updates to `plot_pca()` and gene-level heatmaps.



# bcbioRNASeq 0.0.18 (2017-05-24)

- Simplified count pooling functions.



# bcbioRNASeq 0.0.17 (2017-05-21)

- Reduced number of exports and improved documentation.



# bcbioRNASeq 0.0.16 (2017-05-18)

- Draft migration of [bcbio][] run object into S4 `bcbioRnaDataSet`.
- Created a new variant of `load_run()` that saves to S4 object instead of list.



# bcbioRNASeq 0.0.15 (2017-05-15)

- Reworked and re-organized internal functions.



# bcbioRNASeq 0.0.14 (2017-05-10)

- Defaulted to loading run using project summary YAML file.
- Initial commit of [RMarkdown][] templates (e.g. quality control).
- Added support for dynamic file downloads from [HBC][] website.
- Draft build of website using `pkgdown::build_site()`.



# bcbioRNASeq 0.0.13 (2017-05-08)

- Improved [RDAVIDWebService][] utility functions to work with [dplyr][] 0.6.0.



# bcbioRNASeq 0.0.12 (2017-05-01)

- Reworked metadata and summary metrics functions to obtain information from `project-summary.yaml` saved in the final run directory.



# bcbioRNASeq 0.0.11 (2017-04-27)

- Reduced number of depdencies.
- Initial commit of modified volcano plot from [CHBUtils][] package.
- Internal code updates for upcoming [dplyr][] 0.6.0/[tidyeval][] update.
- Updated [Ensembl][] [biomaRt][] annotations to use live site, currently release 88.



# bcbioRNASeq 0.0.10 (2017-04-19)

- Renamed `import_*` functions to `read_*`.



# bcbioRNASeq 0.0.9 (2017-04-13)

- Consolidated NAMESPACE imports.
- Defaulted to writing count matrices with gzip compression, to save disk space.



# bcbioRNASeq 0.0.8 (2017-04-12)

- Renamed internal parameters for better readability.
- Improved documentation and consolidate functions by group.



# bcbioRNASeq 0.0.7 (2017-04-10)

- NAMESPACE simplification using [basejump][] package.



# bcbioRNASeq 0.0.6 (2017-04-07)

- Reworked handling of plots and tables during knits.



# bcbioRNASeq 0.0.5 (2017-04-06)

- Initial commit of differential expression and gene set enrichment functions.



# bcbioRNASeq 0.0.4 (2017-04-04)

- Added [bcbio][] object integrity checks.
- Improved detection and handling of lane split samples.



# bcbioRNASeq 0.0.3 (2017-03-31)

- Reworked functions to utilize [bcbio][] list object.



# bcbioRNASeq 0.0.2 (2017-03-28)

- Added plotting functions.



# bcbioRNASeq 0.0.1 (2017-03-22)

- Start of package development.
- Initial draft release supporting automatic loading of [bcbio][] run data.



[basejump]: http://steinbaugh.com/basejump
[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[bcbioSingleCell]: https://github.com/hbc/bcbioSingleCell
[bcbioSmallRNA]: https://github.com/lpantano/bcbioSmallRna
[biomaRt]: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
[CHBUtils]: https://github.com/hbc/CHBUtils
[covr]: https://github.com/jimhester/covr
[dplyr]: http://dplyr.tidyverse.org
[Ensembl]: http://www.ensembl.org
[F1000]: https://f1000.com
[HBC]: http://bioinformatics.sph.harvard.edu
[lintr]: https://github.com/jimhester/lintr
[Markdown]: https://daringfireball.net/projects/markdown/syntax
[RDAVIDWebService]: https://bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html
[RMarkdown]: http://rmarkdown.rstudio.com
[testthat]: https://github.com/hadley/testthat
[tidyeval]: http://dplyr.tidyverse.org/articles/programming.html
[Travis-CI]: https://travis-ci.org
[viridis]: https://cran.r-project.org/web/packages/viridis/index.html
