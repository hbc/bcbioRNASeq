# Release notes

## bcbioRNASeq 0.5.1 (UNRELEASED)

Minor changes:

- Updated functional analysis R Markdown template.
- Removed `sessioninfo` as a dependency, switching back to using utils
  `sessionInfo` internally instead.

## bcbioRNASeq 0.5.0 (2022-05-09)

Major changes:

- Starting a new release series requiring R 4.2 / Bioconductor 3.15.
- Minor bug fixes to provide compatibilty with Acid Genomics dependencies.

Minor changes:

- Now using `DataFrame` return from `metrics` by default. Hardened calls to
  `ggplot` using `as.data.frame` coercion.

## bcbioRNASeq 0.4.1 (2022-04-29)

Major changes:

- Now requiring R 4.2 / Bioconductor 3.15.

Minor changes:

- Minor S4 class fixes to provide compatiblity with Bioconductor 3.15.
- Relaxed version pins in `DESCRIPTION` file.

## bcbioRNASeq 0.4.0 (2022-03-11)

Starting a new release series for 2022.

Major changes:

- Package now requires R 4.1+ / Bioconductor 3.14+.
- Package no longer attaches basejump automatically.

Minor changes:

- Removed internal dependencies on `formalsList` parameters defined in basejump,
  which will be removed in a pending release update.

## bcbioRNASeq 0.3.44 (2021-12-14)

Minor changes:

- Reworked the default expected file path of `bcbioRNASeq` object
  (e.g. `"bcb.rds"`), which has been moved back to `data` instead of nested
  inside of `rds` by date (e.g. YYYY-MM-DD).
- Other R Markdown templates, specifically the differential expression and
  functional analysis templates have also been updated to reflect this change.
- Minimum dependencies have been bumped back to R 4.0 / Bioconductor 3.13
  due to current bioconda build issues, which will be resolved in a future
  update of the package.

## bcbioRNASeq 0.3.43 (2021-12-13)

Minor changes:

- Updated dependencies to require Bioconductor 3.14.
- Bug fix for assert check failure of tximport list return. The `infReps`
  slot in list does not return named rows when `txOut` is `TRUE` for
  Bioconductor 3.14 release.

## bcbioRNASeq 0.3.42 (2021-07-27)

Minor changes:

- Updated dependencies, in particular now requiring fixed version of
  DESeqAnalysis that labels genes in `plotVolcano` call as expected.
- Updated comments in functional analysis template.
- Added additional suggested packages that are required for R Markdown
  templates to render when called from a conda environment.

## bcbioRNASeq 0.3.41 (2021-07-21)

Minor changes:

- Improved R Markdown template support, particularly handling of Entrez
  identifier mapping in the functional analysis template. Also updated
  pathview code to timeout on very large pathways, such as general metabolism,
  which can contain hundreds of genes.
- Updated unit tests against F1000 example dataset.
- Removed usage of magrittr pipe in examples, in favor of base R pipe (4.1+).
- Improved internal package name and version handling.
- Harden internal tximport handling code a bit.
- Main `bcbioRNASeq` importer now supports a very minimal dataset with a
  single sample, which isn't common.

## bcbioRNASeq 0.3.40 (2021-02-23)

Maintenance release, incorporating useful changes from the basejump v0.14
release series.

Minor changes:

- Reworked S4 generic imports, migrating some functions to AcidGenerics.
- Reduced the number of dependency packages, since we can offload some of the
  imports to basejump.
- Switched from cli to AcidCLI (via basejump) internally.
- Note that internal GenomicRanges utilities (e.g. `makeGRangesFromEnsembl`)
  need to have `ignoreVersion = TRUE` defined, since this argument changed
  by default in AcidGenomes.

## bcbioRNASeq 0.3.39 (2020-12-03)

Minor changes:

- Bug fix for featureCounts directory structure change in bcbio v1.2.4 release.
- Wrapped DESeq2 normalization calculations (to be defined in `assays` slot)
  in a `tryCatch` call, to more gracefully avoid errors on minimal data sets,
  such as the bcbio-nextgen unit test data.

## bcbioRNASeq 0.3.38 (2020-10-08)

Minor changes:

- Maintenance release, updating minimum dependency package versions.

## bcbioRNASeq 0.3.37 (2020-09-15)

Minor changes:

- Hardened plotting functions against bcbioRNASeq objects created with
  `fast = TRUE`.
- `show`: Added display of `fast: TRUE` when applicable.
- Documentation updates.

## bcbioRNASeq 0.3.36 (2020-09-11)

Minor changes:

- Decreased R package dependency versions a bit, so we can update the bioconda
  recipe successfully. Some Bioconductor packages don't have point release
  updates on bioconda.

## bcbioRNASeq 0.3.35 (2020-07-24)

Minor changes:

- Maintenance release, updating minimum R dependency to 4.0.

## bcbioRNASeq 0.3.34 (2020-07-08)

Minor changes:

- Temporarily disabled `plotPCACovariates`, which causing CI checks to error due
  to a dplyr bug in DEGreport. Will reenable when this gets fixed and is updated
  on Bioconductor.

## bcbioRNASeq 0.3.33 (2020-05-11)

Minor changes:

- `bcbioRNASeq` generator now has improved message support for selected samples.
  Internally, this is managed by the cli package.
- `updateObject`: Now supporting `verbose` argument, which defaults to `FALSE`.

## bcbioRNASeq 0.3.32 (2020-02-24)

Minor changes:

- `bcbioRNASeq`: Removed now defunct `spikeNames` argument. Refer to
  `makeSummarizedExperiment` documentation and release notes in basejump
  package for details.

## bcbioRNASeq 0.3.31 (2020-02-19)

Minor changes:

- Relaxed validity checks against `dataVersions` and `programsVersions` in
  `metadata` slot, for compatibility with `updateObject` upgrade path from
  v0.2.9 objects, which are currently generated by bcbio-nextgen. Thanks
  @sjmonkley for catching this. See issue #146.
- Changed license from MIT to GPL-3.
- Updated documentation.

## bcbioRNASeq 0.3.30 (2020-01-20)

Major changes:

- Improved internal handling of `bcbioRNASeq` to `DGEList` coercion, for handoff
  to edgeR and limma-voom. Now the code handles `countsFromAbundance` settings
  more intelligently and informs the user whether downstream edgeR and/or
  limma-voom handoff is appropriate.
- Now exporting `as.DESeqDataSet`, `as.DESeqTransform`, and `as.DGEList` S4
  generic methods, which support arguments. Note that quick coercion using
  only `as` (see `coerce` documentation from Bioconductor) does not support
  arguments.

Minor changes:

- Fixed support for bcbio fastrnaseq pipeline. Thanks to @kokyriakidis for
  pointing this out. Added code coverage support for an example run produced
  via this pipeline.
- Updated basejump dependencies, namely renaming bioverbs to AcidGenerics.
- Now using cli package for improved messages.

## bcbioRNASeq 0.3.29 (2019-10-30)

Minor changes:

- Added support for Bioconductor 3.10.
- Updated internal `DataFrame` validity checks to also check for `DFrame` class,
  which has changed in latest S4Vectors 0.23 release.
- Updated basejump dependency versions.

## bcbioRNASeq 0.3.28 (2019-09-17)

Minor changes:

- Reworked quality control plots to use "labels" instead of "title".
- Updated R Markdown templates and added render code coverage in testthat-extra.
- Updated vignette to use shared bibliography.bib file from basejump.

## bcbioRNASeq 0.3.27 (2019-09-09)

Minor changes:

- `bcbioRNASeq`: Internal generator now calls `importSampleData` using
  `pipeline = "bcbio"` argument, so we don't run into breaking changes when
  handling user metadata in a future basejump update.
- Updated basejump dependency versions.

## bcbioRNASeq 0.3.26 (2019-08-27)

Minor changes:

- Now requiring R 3.6, as we near a candidate Bioconductor release.
- Updated basejump dependencies.

## bcbioRNASeq 0.3.25 (2019-08-20)

Major changes:

- Made subsetting factor for extract (`[`) operations.
- `updateObject`: Improved `rowRanges` handling.

Minor changes:

- `bcbioRNASeq` generator now captures call using `standardizeCall` instead of
  `match.call`, and automatically expands all default arguments.
- Improved NAMESPACE and basejump dependencies.
- Made vsn and DEGreport packages optional, to reduce the number of dependencies
  required for installation.
- `plotPCACovariates`: Reworked internal code to not depend on dplyr.
- Removed internal code dependencies on magrittr pipe.

## bcbioRNASeq 0.3.24 (2019-08-12)

Minor changes:

- Simplified internal object creation in `bcbioRNASeq`, `[`, and `updateObject`.
- Improved message consistency.

## bcbioRNASeq 0.3.23 (2019-07-29)

Major changes:

- `updateObject`: Worked internal approach to be compatible with upcoming
  Bioconductor 3.10 release. Has to rethink the strategy to handling assays
  internally, which can no longer be easily coerced from assays slot in the
  `bcbioRNASeq` object, due to changes in `ShallowSimpleListAssays` handling.
  This should return an object in exactly the same manner as previous releases.

## bcbioRNASeq 0.3.22 (2019-07-23)

Major changes:

- `bcbioRNASeq`: Now supporting `fast` argument, which enables the user to skip
  any internal DESeq2 calculations. This will skip generation of `normalized`,
  `vst`, and `fpkm` matrices. Note that automatic `rlog` matrix generation is
  no longer supported, since it's often too slow, especially for large datasets.
  If you want to generate rlog counts, simply coerce the bcbioRNASeq object
  to a DESeqDataSet, then run the `rlog` function.
- Aligned counts are now saved inside the object when a pseudoaligner (e.g.
  salmon, kallisto) are used to define the primary `counts` assay. These values
  are saved in an `aligned` matrix.

New functions:

- `pseudoVsAligned`: Visually inspect the correlation of pseudoaligned counts
  to the aligned counts. Note that default bcbio genomes use legacy UCSC builds,
  where as the pseudoaligners generate values from newer Ensembl transcript
  annotations.

Minor changes:

- Updated package dependency versions.

## bcbioRNASeq 0.3.21 (2019-07-17)

Major changes:

- Initial support for `bcbioRNASeq` to `DGEList` (edgeR) object coercion.
- `relativeLogExpression`: Improved internal code, using edgeR
  `calcNormFactors` approach.

Minor changes:

- Updated basejump dependency to v0.10.11.
- Improved Travis CI docker config to check coverage and lints.
- Improved Bioconductor installation instructions.

## bcbioRNASeq 0.3.20 (2019-05-29)

Minor changes:

- Added support for size factor adjusted counts as "sf" in `counts` call.
- Improved internal trans axis labeling handling for DESeqTransform counts.
- `plotMeanSD`: Improved default colors. Improved internal non-zero count
  handling. Reworked plot headers.

## bcbioRNASeq 0.3.19 (2019-05-05)

Major changes:

- Now pinned to R >= 3.5.

Minor changes:

- Improved Travis CI and AppVeyor CI configuration.

## bcbioRNASeq 0.3.18 (2019-04-25)

Minor changes:

- S4 generic reexport documentation fixes.

## bcbioRNASeq 0.3.17 (2019-04-22)

Minor changes:

- Bar plots now show the first sample at the top of the Y axis when
  `flip = TRUE`. This is more human readable.
- Quality control template now uses `bcbiornaseq_file` as main input param.
- Bug fix: added back internal `.normalizedTrans` function required for some
  quality control plots.

## bcbioRNASeq 0.3.16 (2019-04-17)

Minor changes:

- Tweaks to quality control template.
- Switched Travis CI configuration to use `rnaseq` Docker image.

## bcbioRNASeq 0.3.15 (2019-04-11)

Minor changes:

- Ensuring that all [basejump][] functions necessary for the quality control
  template get reexported.
- Renamed R Markdown template files and `travis-render.sh` script into kebab
  case, rather than snake case.
- Documentation improvements (using `inheritParams`).

## bcbioRNASeq 0.3.14 (2019-04-01)

Moved forked development repository to [Acid Genomics][].

Major changes:

- Switched unit testing data URL from seq.cloud to acidgenomics.com.

Minor changes:

- Reworked internal reexports for assignment generics. See `counts` for example.
- Miscellaneous documentation improvements.
- Switched from `basejump` prefix to `acid`, where applicable.
- Reexported essential functions from basejump required for QC template to
  render without basejump being attached.

## bcbioRNASeq 0.3.13 (2019-03-29)

Minor changes:

- Added additional checks to `bcbioRNASeq` validity method.
- Draft update to quality control template, to support parameterized
  `bcbioRNASeq` object input. This is likely not necessary and will be removed
  in a future update.
- Added additional checks for `bcbioRNASeq` to `DESeqDataSet` coercion.
- Improved internal `perMillion` argument code handling inside plotting
  functions: `plotMappedReads`, `plotTotalReads`.
- Removed `reexports.R`. No longer reexporting `Gene2Symbol` or magrittr pipe.

## bcbioRNASeq 0.3.12 (2019-03-18)

Major changes:

- `sampleData`: Moved method support to [basejump][] package.

Minor changes:

- Resaved example data.
- Sample metadata chunk in quality control template now returns `DataFrame`.
- Improved some deprecation messages.
- CI configuration improvements.

## bcbioRNASeq 0.3.11 (2019-02-12)

Minor changes:

- Renamed `plotGene` to `plotCounts`, to match consistency in [basejump][] and
  [DESeqAnalysis][] package updates.
- Updated dependencies to provide improved backwards compatibility support for
  R 3.4 / Bioconductor 3.6.

## bcbioRNASeq 0.3.10 (2019-02-01)

Minor changes:

- Travis CI build configuration and documentation fixes.

## bcbioRNASeq 0.3.9 (2019-01-23)

Minor changes:

- Bug fixes for `interestingGroups` handling in plotting functions. Needed
  to add back an `interestingGroups<-` assignment call to slot the object.
- Removed `sanitizeRowData` and `sanitizeSampleData` from imports.

## bcbioRNASeq 0.3.8 (2019-01-17)

Minor changes:

- `.tximport`: Use `makeNames` to check for valid names internally.
- Improved documentation.
- Travis CI configuration fixes to get build checks to pass.

## bcbioRNASeq 0.3.7 (2019-01-13)

Major changes:

- Overhauled validity method for `bcbioRNASeq`. Improved check steps using `ok`
  method, similar to approach in goalie.
- Added back `plotPCACovariates` method support.

Minor changes:

- Updated pkgdown configuration. Removed `alphaSummary`, `contrastName`,
  `export`, `markdown`, `sampleData`, `topTables`
- goalie package has been updated to support `nullOK` mode where applicable.
  This is particularly use for formals that use `NULL` by default. See
  `bcbioRNASeq` generator changes, for example.
- `prepareRNASeqTemplate` example is causing build checks to fail, so disable.
- Improved CI configuration.

## bcbioRNASeq 0.3.6 (2018-12-12)

This is the first release that switched over to using goalie for assert checks
instead of assertive.

Major changes:

- `bcbioRNASeq` generator has been overhauled to use `assert` for assert checks.

Minor changes:

- Ensure title is quoted in R Markdown.
- Switched to using `validate` in place of `validate_that` in validity checks.
- Miscellaneous documentation improvements.
- Now using `match.arg` internally to check validity of `countsFromAbundance`
  argument passthrough to `tximport`.

## bcbioRNASeq 0.3.5 (2018-12-01)

Major changes:

- Added additional improvements to quality control template. Ensuring that setup
  chunk never caches. bcbioRNASeq object should be called using
  `params$object_file`. This convention will be used across all other bcbio R
  packages. Reworked and simplified `plotCountsPerGene` section of QC template.

Minor changes:

- Added early return method for objects without modification extract method.
  Might want to reconsider this approach if the user just wants to reorder
  samples. We're checking for unmodified return by looking for identical
  raw counts matrix before and after extraction call.
- `plotCountsPerBiotype`: Ensure `normalized` argument uses "tpm" first.
- `plotDispEsts`: Switch to using `hasUniqueCols` internally instead of
  `areSamplesUnique`.

## bcbioRNASeq 0.3.4 (2018-11-29)

Major changes:

- Improved and simplified quality control template. Removed Dropbox mode (which
  may be added back in a future update). Identical samples detection is now
  supported automatically in the plotting functions, which was previously
  defined in the QC template.
- Improved subset sample handling inside extraction method. Extraction method
  now imports useful `relevelRowRanges` and `relevelColData` utilities from
  basejump.

Minor changes:

- Removed `sampleData<-` reexport.
- Documentation fixes and improvements, particularly to the `bcbioRNASeq`
  generator function.
- Reworked internal gffFile handling inside `bcbioRNASeq` generator function.
- Simplified passthrough in internal `.new.bcbioRNASeq` call.
- `plot5Prime3PrimeBias` now uses `color` argument instead of `fill`.
- `plotDispEsts` now checks for duplicate samples and early returns.
- Consistently use British spelling variants (colour instead of color) for
  ggplot2 functions. Our formals will stick with American English though.
- `plotTotalReads`: Bumped recommended default to 20 million from 10 million.

## bcbioRNASeq 0.3.3 (2018-11-25)

Minor changes:

- Added additional library loads to `_setup.R`: rmarkdown, basejump, DESeq2.
- `plotGene`: Updated working example.
- Improved reexport of `sampleData<-` generic.
- Updated Travis CI configuration.

## bcbioRNASeq 0.3.2 (2018-11-19)

Major changes:

- Splitting out new `DESeqAnalysis` S4 class to a separate package.

Deprecated functions and methods:

- Functions removed from `NAMESPACE`: magrittr pipe (`%>%`), `DESeqAnalysis`,
  `DESeqResultsTables`, `alphaSummary`, `contrastName`, `export`, `markdown`,
  `plotDEGHeatmap`, `plotDEGPCA`, `plotMA`, `plotMeanAverage`, `plotVolcano`,
  `topTables`.
- `alphaSummary`: Migrated `DESeqDataSet` method to [DESeqAnalysis][] package.

Example data:

- Deleted `deseq` example, which has migrated to [DESeqAnalysis][] package.

Minor changes:

- Reorganized internal assert checks.
- Miscellaneous documentation improvements.
- Improved consistency of global formals (using `getOption`) via `formalsList`
  global variable.

## bcbioRNASeq 0.3.1 (2018-11-14)

Major changes:

- Improved validity method for `bcbioRNASeq` S4 class.
- `extract` (`[`) method: Improved internal handling of `rowRanges` and
  `colData`. Now using `I` internally for complex S4 columns in `mcols` for
  `rowRanges`.
- Reworked imports back to basejump instead of attempting to use subpackages:
  basejump.annotations, basejump.assertions, basejump.classes, basejump.coerion,
  basejump.developer, basejump.experiment, basejump.generics, basejump.globals,
  basejump.io, basejump.markdown, basejump.plots, basejump.sanitization. Rethink
  this approach in a future update.

Minor changes:

- Resaved example `bcb` dataset.
- Miscellaneous internal assert check fixes and updates.
- Miscellaneous documentation improvements.
- Added internal `Rle` global variable. May want to avoid this approach because
  `Rle` is a defined function in S4Vectors, however.
- Offloaded some global documentation params to basejump.
- Reworked some global params set by `getOption` in formals.

## bcbioRNASeq 0.3.0 (2018-11-06)

Development fork of `hbc/bcbioRNASeq`. Code development will proceed here on the
fork until a stable release is ready to merge back at official HBC repo, to
avoid any disruption for active RNA-seq experiments.

New S4 classes:

- Added new `DESeqAnalysis` S4 class, which helps containerize up a DESeq2
  analysis, often requiring corresponding `DESeqDataSet`, `DESeqTransform`, and
  `DESeqResults` objects. Note that this class may be defined in a separate
  package in a future update.

Code migration to basejump:

- `aggregateReplicates` has moved to [basejump][], and now works primarily on
  `SummarizedExperiment` method, instead of `bcbioRNASeq`.
- Offloaded a number of S4 generics to basejump infrastructure, which are now
  defined in the bioverbs package: `alphaSummary`, `contrastName`,
  `plot5Prime3PrimeBias`, `plotCountDensity`, `plotCountsPerGene`,
  `plotDEGHeatmap`, `plotDEGPCA`, `plotExonicMappingRate`, `plotGenderMarkers`,
  `plotGeneSaturation`, `plotGenesDetected`, `plotIntronicMappingRate`.
- `metrics` method, which returns a `tbl_df` of `colData`, has been moved to
  basejump.
- `plotCorrelationHeatmap` code is now primarily defined in basejump as a
  `SummarizedExperiment` method. bcbioRNASeq now calls this method internally
  but keeps support for `normalized` counts argument.
- Removed `tpm` method support, which is now handled by `SummarizedExperiment`
  method support in basejump.

S4 method reworks:

- Reworked internal code, splitting out S4 methods into internal functions
  (e.g. `counts.bcbioRNASeq`). This step is recommended by [Bioconductor][], and
  can make troubleshooting easier.
- `export`: Reworked S4 methods for `DESeqResults` and `DESeqResultsTables` S4
  classes.
- `extract` (i.e. `[`) method: Switched `transform` formal to `recalculate`.
  Added additional assert checks. Operation no longer changes the saved
  `version` of bcbioRNASeq, slotted into `metadata`.
- `interestingGroups` formal is now set `NULL` instead of missing for plotting
  functions.
- `markdown`: Initial support of new S4 generic, which returns a markdown table
  of sample metadata.
- `plot5Prime3PrimeBias`: Added an addition grep matching step to handle
  differential variations on the `x5x3Bias` column in `colData`.
- `plotCountDensity` S4 method has been removed. Use `plotGeneSaturation`
  generic with `geom = "density"` argument instead.
- `plotCountsPerBiotype`: New S4 method support, which calls basejump
  `SummarizedExperiment` method.
- `plotDEGHeatmap`: Draft update of `DESeqResults` method.
- `plotDEGPCA`: Draft update of `DESeqResults` method.
- `plotDispEsts`: Split out and reworked internal method definition.
- `plotGenesDetected`: This code has been moved to basejump and is now defined
  against `SummarizedExperiment`.
- `plotMA`: Reworked internal `DESeqResults` method code. Added draft support
  for new `DESeqAnalysis` S4 class. Keeping `plotMeanAverage` reexported, since
  this function is used in the F1000 workflow paper.
- `plotMappingRate`: Set the default limit back down to 0.7.
- `plotMeanSD`: Reworked internal `DESeqDataSet` code. Added additional assert
  checks and changed the default handling of `vst` and `rlog` internal formals.
- `plotPCA`: Draft update of `SummarizedExperiment`. Will offload this code to
  basejump in a future update. Switched `return` argument to support `tibble`
  and S4 `DataFrame` -- note change from previous support for `data.frame`.
- Deleted `plotPCACovariates` method, since DEGreport is failing. Will add back
  support for this generic in a future update.
- `plotQC`: Initial method support for `bcbioRNASeq`.
- `plotTotalReads`: Reworked and tightened up internal ggplot2 code.
- `plotVolcano`: Draft update to `DESeqResults` method. Will rework this code
  with a `DESeqAnalysis` approach in a future update.
- `relativeLogExpression`: Added S4 method support.
- `resultsTables`: Removed this approach to subsetting `DESeqResults`. This will
  be a rework in the upcoming [DESeqAnalysis][] package.
- `sampleData`: Simplified internal code, which takes advantage of updates to
  `SummarizedExperiment` defined in basejump.
- `updateObject`: Reworked internal metadata handling.

Example data reorganization:

- Renamed and simplified example datasets. Now only exporting `bcb`
  (`bcbioRNASeq`) and `deseq` (`DESeqDataSet`) examples.
- `bcb_small`, `dds_small`, `gender_markers`, and `res_small` have been removed.
- Updated and consolidated corresponding `data-raw/` scripts.
- Reworked internal `extdata/bcbio/` example datasets. Testing out both
  `2017-05-23_rnaseq` and `2018-03-18_GSE67267-merged` data.

R Markdown templates:

- Stripped down supported R Markdown templates. Currently simplified code base
  to support Quality Control markdown (`quality_control`). Differential
  expression and functional analysis templates have been temporarily removed
  but will be re-added in the future.

Minor changes:

- Reorganized S4 class definitions and validity checks in `AllClasses.R`.
  Previously, this code was split out per S4 (e.g. `bcbioRNASeq-validity.R`).
- Moved S4 coercion methods from `coerce-methods.R` to `as-methods.R`.
- Added internal assert checks in `assert-internal.R`.
- Consolidated internal DESeq2 code into `DESeq2-internal.R`.
- Consolidated and defined additional global variables into `globals.R`.
- Reorganized internal `.meltCounts` code.
- Consolidated global package parameter arguments into `params.R`.
- `bcbio` prefix has been renamed to `basejump`, where applicable. Refer to
  global function params defined with `getOption`, and internal ggplot2 utility
  functions, such as `basejump_geom_abline`.
- Draft support for `countsFromAbundance` selection in internal `tximport` call.
  This is useful for some transcript-level analyses.

Deprecations:

- Removed defunct `download` function.
- New deprecations: `plotCountDensity`, `plotPCACovariates`, `resultsTables`.

## bcbioRNASeq 0.2.10 (2019-05-29)

Minor changes:

- Update to support basejump v0.10+ release series.

## bcbioRNASeq 0.2.9 (2019-05-03)

Major changes:

- Fixed a bug in internal tximport code (refer to `tximport-internal.R` file)
  where sample names can get associated with the wrong samples in experiments,
  causing samples with numbers roll over from 1 -> 10. Thanks to Richard from
  AstraZeneca for catching this.

## bcbioRNASeq 0.2.8 (2019-01-28)

- Handle correlation of PCA covariates when a numerical covariate has no
  variation.
- Handle `NA` in ribosomal RNA calculations.

## bcbioRNASeq 0.2.7 (2018-08-22)

Minor changes:

- Bug fix for `as` coercion method support. Need to ensure
  `exportMethods(coerce)` is included in `NAMESPACE` file, otherwise
  `bcbioRNASeq` to `DESeqDataSet` coercion using `as(bcb, "DESeqDataSet")` won't
  work when called from an Rscript without the package library loaded. Thanks
  @roryk for noticing this.
- Only run rRNA plots if rRNA is detected.

## bcbioRNASeq 0.2.6 (2018-08-19)

Major changes:

- Added `transform = TRUE` argument to `[` extraction method, allowing the
  user to skip automatic [DESeq2][] transformations, which can be CPU intensive
  for large datasets.
- Switched back to using `plotMA` instead of `plotMeanAverage`. An MA-plot
  by definition is not a "Mean Average" plot, so this function name is
  misleading. We will keep the `plotMeanAverage` working but it is now
  soft deprecated.
- `plotGeneSaturation` now supports `label` argument, similar to `plotPCA`.

Minor changes:

- Improved internal handling during `tximport` call to handle transcript
  version mismatch with tx2gene data.frame. This can result if the bcbio
  pipeline is using an old genome build.
- Ensure `genomeBuild` is detected from [AnnotationHub][] `rowRangesMetadata` if
  applicable, and not left `NULL` in the metadata.
- Updated internal code to use `aes` instead of `aes_string`, which uses
  tidyeval and quasiquotation.
- `plotGene`: reduced the number of `return` options to simply "facet" and
  "wide". Previously, this also supported "grid", "list", and "markdown", but
  these were removed because they are not frequently used.
- `plotGene`: Switched back to internal `lapply` call instead of using
  `BiocParallel::bplapply`. This doesn't always work perfect in an HPC
  environment (e.g. HMS O2 cluster).
- Soft deprecated `plotMeanAverage` in favor of `plotMA`.
- All S4 generics with method support are exported for easier lookup in the
  reference documentation.
- All function reexports have been consolidated in [bcbioBase][] package.

## bcbioRNASeq 0.2.5 (2018-06-21)

Minor changes:

- Enable support [bcbio][] integration, by modifying `bcbioRNASeq` constructor
  to work with minimal bcbio test data.
- Switched from internal usage of `aes_` in favor of consistent usage of
  `aes_string`. This will make the transition to [ggplot2][] v2.3.0 easier
  in a future update.

## bcbioRNASeq 0.2.4 (2018-05-24)

Major changes:

- `aggregateReplicates` support has been added back. This function returns
  a `RangedSummarizedExperiment` instead of a `bcbioRNASeq` object, containing
  only an aggregate raw counts matrix in the `counts` slot of `assays`.
- The functional analysis [R Markdown][] template has been reworked to use
  `dds_file` and `organism` as new parameter arguments. We've reduced the number
  of parameters required here to run [clusterProfiler][].
- Made `alphaSummary` defunct for `bcbioRNASeq` object, in favor of
  `DESeqDataSet` only. This function is only useful when a proper design formula
  has been defined.

Minor changes:

- `metrics` now contains an informative error for datasets that were analyzed
  using the `fast-rnaseq` [bcbio][] pipeline.
- `DESeqDataSet` coercion from `bcbioRNASeq` object doesn't attempt to run
  `DESeq` command any more, which was unnecessary and improves speed.
- `bcbioSingleCell` constructor now supports `censorSamples` parameter. This
  is useful for removing known poor quality samples upon loading.
- [ggplot2][] color and fill palettes are now set `NULL` in the quality control
  functions. This behavior doesn't change the appearance of the plot colors,
  which will still default to `ggplot2::scale_colour_hue` or
  `ggplot2::scale_fill_hue`. The upcoming [ggplot2][] v2.3.0 update supports
  global options for color and fill palettes, so these parameters may be
  deprecated in a future release.
- Reworked the internal code for `topTables`.

Infrastructure changes:

- Added macOS testing to Travis CI build checks.
- Fixed [clusterProfiler][] compilation error on Travis CI by installing
  `libudunits2-dev` (Linux).

## bcbioRNASeq 0.2.3 (2018-05-10)

Major changes:

- Now recommending variance stabilizing transformation (`vst`) over `rlog`
  counts by default in plots, where applicable.

Minor changes:

- Tweaked Rory's biotype plots in QC report to match formatting conventions in
  the package. These plots are now colored.
- Added `plotDEGPCA` to default differential expression [R Markdown][]
  template.
- `colData` factors are correctly releveled upon object subset with `[`. This
  helps avoid unwanted downstream errors when creating a `DESeqDataSet` and
  running differential expression with [DESeq2][].
- Recommending `facet` return method by default for `plotGene`. Updated the
  working example to reflect this.
- `metrics` now returns `interestingGroups` column.
- `sample` label has been removed from axis title for QC plot functions.
- Now using shared [ggplot2][] convneience functions from [bcbioBase][] 0.2.10:
  `bcbio_geom_abline`, `bcbio_geom_label`, and `bcbio_geom_label_repel`.
  These are also used by [bcbioSingleCell][] for improved graphical consistency.
- Removed unused internal legacy [ggplot2][] code.
- Increased [DEGreport][], [DESeq2][], and [tximport][] dependency requirements.

## bcbioRNASeq 0.2.2 (2018-04-26)

Minor changes:

- Split out assertive imports so we can pin on [bioconda][].
- Improved package documentation.
- Improved label consistency in `plotPCA` functions to match
  `plotMeanAverage` and `plotVolcano`.
- Improved automatic title labeling in `plotDEGPCA`, matching the other DEG
  functions. Also added directionality to `plotDEGPCA`.
- Added `DESeqDataSet` method support to `plotCorrelationHeatmap`, using the
  normalized counts.
- `reusltsTables` now writes local files to `tempdir` when [Dropbox][] mode
  is enabled using `dropboxDir`.

## bcbioRNASeq 0.2.1 (2018-04-24)

Last set of code fixes before F1000v2 resubmission.

Major changes:

- Added `rle` return support for `counts`, which are calculated on the fly.
- Added `transgeneNames` and `spikeNames` support to `loadRNASeq` function.
- `loadRNASeq` now supports `organism = NULL` again, for datasets with poorly
  annotated genomes.
- Primary `assay` containing raw counts is now named `counts` instead of
  `raw`, for consistency with other `SummarizedExperiment` objects (e.g.
  `DESeqDataSet`) and the [bcbioSingleCell][] S4 class definition.
- Improved internal code for `plotGene` and `plotGenderMarkers`.

Minor changes:

- Improved [AppVeyor CI][] support to test against bioc-devel using [R][] 3.5.
- Improved support and unit testing for `updateObject` method.
- [DESeq2][] normalized counts are always slotted in `assays`, even when rlog
  and vst transformations are skipped.
- Exporting `[[<-`, `assays<-`, `colData<-`, `interestingGroups<-`, and
  `metadata<-` assignment methods, to avoid unwanted coercion to
  `SummarizedExperiment`. Objects extending `RangedSummarizedExperiment`
  shouldn't be doing this, so we may need to file a bug report with Bioconductor
  or check our class definition in the package.
- Now importing specific functions from S4Vectors and methods rather than
  importing everything.
- Switched back to using `stop`, `warning` and `message` rather than the
  alternate [rlang][] functions `abort`, `warn`, and `inform`.
- Objects with invalid metadata now print which slots are invalid to the
  console.

## bcbioRNASeq 0.2.0 (2018-03-22)

Major changes:

- `bcbioRNASeq` S4 class object is now extending `RangedSummarizedExperiment`
  instead of `SummarizedExperiment`. Consequently, the row annotations are now
  stored in the `rowRanges` slot as `GRanges` class, instead of in the
  `rowData` slot as a `DataFrame`. The `rowData` accessor still works and
  returns a data frame of gene/transcript annotations, but these are now
  coerced from the internally stored `GRanges`. The `GRanges` object is
  acquired automatically from [Ensembl][] using `basejump::ensembl`. By
  default, `GRanges` are acquired from [Ensembl][] using [AnnotationHub][] and
  [ensembldb][]. Legacy GRCh37 genome build is supported using the
  [EnsDb.Hsapiens.v75][] package.
- `assays` now only slot matrices. We've moved the [tximport][] data from
  the now defunct `bcbio` slot to assays. This includes the `lengths` matrix
  from [tximport][]. Additionally, we are optionally slotting [DESeq2][]
  variance-stabilized counts ("`rlog`", `"vst"`). [DESeq2][] normalized counts
  and [edgeR][] TMM counts are calculated on the fly and no longer stored
  inside the `bcbioRNASeq` object.
- `colData` now defaults to returning as `data.frame` instead of
  `DataFrame`, for easy piping to [tidyverse][] functions.
- `bcbio` slot is now defunct.
- FASTA spike-ins (e.g. EGFP, ERCCs) can be defined using the `isSpike`
  argument during the `loadRNASeq` data import step.
- Melted counts are now scaled to log2 in the relevant quality control
  functions rather than using log10. This applies to `plotCountsPerGene` and
  `plotCountDensity`. Note that we are subsetting the nonzero genes as
  defined by the raw counts here.
- Simplified internal `tximport` code to no longer attempt to strip
  transcript versions. This is required for working with _C. elegans_
  transcripts.
- Minimal working example dataset is now derived from GSE65267, which is also
  used in the F1000 paper.
- Added `as(object, "DESeqDataSet")` coercion method support for
  `bcbioRNASeq` class. This helps us set up the differential expression
  analysis easily.
- `counts` function now returns DESeq2 normalized counts (`normalized = TRUE`)
  and edgeR TMM counts (`normalized = "tmm"`) on the fly, as suggested by the
  F1000 reviewers.
- Design formula can no longer be slotted into `bcbioRNASeq` object, since
  we're not stashing a `DESeqDataSet` any more.
- Updated Functional Analysis R Markdown template.

Minor changes:

- `validObject` is now required for all plotting functions. This check is
  also called in the R Markdown template. Legacy objects can be updated using
  `updateObject`.
- `metrics` now returns columns sorted alphabetically.
- Added `contrastName` as a generic function.
- `plotDEGHeatmap` and `plotDEGPCA` generics no longer have `counts`
  defined in the signature. The `counts` argument is now only defined in the
  methods.
- `prepareRNASeqTemplate` has been converted from a generic to a standard
  function.
- Improved `metadata` validity checks.
- `plotCorrelationHeatmap` matrix method has been moved to [basejump][]
  package, for improved consistency with the other heatmap code.
- `plotGenderMarkers` internal code has been reworked to match
  `plotGene`.
- Default `plotMA` appearance has changed, providing a line at the 0
  y-intercept, similar to `DESeqDataSet` method.
- Internal example datasets have been renamed (e.g. `bcb_small` instead of
  `bcb`).
- Added AppVeyor CI support for code testing on Windows.
- Made Travis CI checks stricter, added `BiocCheck`.
- Internal `.sampleDirs` code is now exported in [bcbioBase][] as a
  generic.
- `gene2symbol` and `interestingGroups` method support are now defined
  for `SummarizedExperiment` in the [bcbioBase][] package.

Updating legacy objects (pre v0.2.0):

- Use `updateObject` in combination with the `rowRanges` argument, which
  requires a `GRanges` object. `GRanges` can be obtained from [Ensembl][] using
  the `basejump::ensembl` function or the [ensembldb][] package.

Deprecations:

- `bcbio` slot is now defunct, since we have moved all data into the
  `SummarizedExperiment` container.
- Deprecated `plot5x3Bias` in favor of `plot5Prime3PrimeBias`. This is
  less confusing as to what this function plots.
- `flatFiles` has been deprecated in favor of `as(object, "list")` coercion
  method. See [bcbioBase][] package for `SummarizedExperiment` method
  support.
- Defunct: `design`, `download`, `meltLog10`, `txi`.
- Legacy `bcbioRNADataSet` method support has been removed.

## bcbioRNASeq 0.1.8 (2018-04-03)

- Bug fix for `gene2symbol` argument not renaming rows in `plotDEGHeatmap`.

## bcbioRNASeq 0.1.7 (2018-02-28)

- Bug fix for `[` subset method dropping metrics in metadata.
- Simplified unit testing for Dropbox mode enabled in `resultsTables`.

## bcbioRNASeq 0.1.6 (2018-02-20)

- Bug fix for gene-to-symbol mappings in `plotDEGHeatmap`.
- Added support for quickly plotting differentially expressed genes (DEG) in
  a PCA plot with `plotDEGPCA`.
- Added support for Dropbox shared links to `resultsTables`, for use with
  the Stem Cell Commons database.
- Added assert checks internally for all functions.
- Improved internal code for `plotGene` and `plotGenderMarkers` to run
  faster.
- Deprecated data frame methods based on metrics for QC functions.

## bcbioRNASeq 0.1.5 (2018-01-31)

- Import shared dependency functions from bcbioBase instead of basejump.
- Added method support for `selectSamples`.
- `organism` and `genomeBuild` parameters are now user-definable in the main
  `loadRNASeq` import function.
- Fixed gene subsetting method on S4 object, which handles genes using
  `intersect` in the featureCounts matrix.
- Removed internal `aggregateReplicates` code. This needs to be reworked
  and added back in a future release.
- Improve method for handling a missing normalized counts matrix in the
  assays slot. This can occur when the user opts to skip the CPU-intensive
  DESeq2 normalizations.
- Improved internal code for the quality control functions. Improved the `if`
  statements to be more class specific.
- Renamed `plotCorrelationHeatmap` `transform` argument to `normalized`,
  for consistency with the `counts` generic.
- Added `title` support to plots, where applicable.
- Updated internal code for `plotDEGHeatmap`.
- Updated internal marker handling code for `plotGenderMarkers`.
- `resulsTables` function now defaults to `summary = TRUE`.

## bcbioRNASeq 0.1.4 (2018-11-27)

- Migrated all basejump function imports to bcbioBase package.

## bcbioRNASeq 0.1.3 (2017-12-03)

- Combined examples (`bcb`, `dds`, `res`, etc.) into a single `examples`
  object. This helps avoid accidental use of example `bcb` in an analysis.
- Moved ggplot imports from `internal-ggplot.R` to above each function.
- Renamed `maxSamples` parameter in `loadRNASeq` to `transformationLimit`.
  If there are more samples than this limit, then the DESeq2 transformations
  will be skipped. In this case, `rlog` and `vst` will not be slotted into
  `assays`.
- Added a colData sanitization step in `loadRNASeq` to ensure rows are in
  the same order as the columns in the counts matrix. Otherwise, DESeq will
  report an error at the `DESeqDataSetFromTximport` step. We're also ensuring
  the factor levels get updated here.
- Now using `glimpse` instead of `str` in examples, where applicable.
- Added `colData<-` assignment method support. This requires a `DataFrame`
  class object. Upon assignment, the internal colData at
  `bcbio(object, "DESeqDataSet")`, `assays(object)[["rlog"]]` and
  `assays(object)[["vst"]]` are also updated to match.
- Initial assignment support for `design`, which will update the internal
  DESeqDataSet.
- Added method support for `gene2symbol` generic, which will now return a 2
  column `data.frame` with `ensgene` and `symbol` columns. This is helpful for
  downstream gene to symbol mapping operations.
- Added working example for `interestingGroups<-` in the documentation.
- Added some code to improve factor releveling, where applicable. See
  `internal-meltLog10.R` for example.
- Now explicitly defining the custom color palettes (e.g.
  `viridis::scale_fill_viridis(discrete = TRUE)`. This makes it clearer to the
  user in the documentation where these palettes are located.
- Improved axis label support in `plotGene`.
- `plotHeatmap` now uses internal `gene2symbol` mappings from stashed
  annotable, instead of always querying Ensembl. The user can define custom
  mappings with the `gene2symbol` argument, if desired.
- `plotPCA` now supports custom color palettes. The `shapes` parameter has
  been removed because it doesn't work well and is limited to datasets with few
  samples. This behavior matches the PCA functionality in DESeq2.
- Improved internal code for `plotVolcano`. Added support for `gene2symbol`
  argument, like in `plotHeatmap`. If left missing, the function will query
  Ensembl for the gene2symbol mappings. We're now using `data` instead of
  `stats` as the main data source.
- Improved legibility of subset method code.
- Added some additional reexports, which are used for the package
  documentation and website.
- Simplified legacy object coercion method code.
- Updated Bioconductor installation method code. We're now using the
  `dependencies` argument, which allows for automatic install of suggested
  packages along with imports.

## bcbioRNASeq 0.1.2 (2017-11-08)

- Updated package imports to match Bioconductor 3.6.
- Added support for interesting groups assignment with `interestingGroups<-`.
- Renamed `plotGeneHeatmap` to simply `plotHeatmap`.
- Added gender marker support for _Homo sapiens_.
- Improved support for multiple interesting groups in quality control plots.
  Now interestingGroups is defined as a column in the metrics data.frame that
  is used to specify the plot color/fill. This matches the convention in the
  bcbioSingleCell 0.0.22 update.
- Sample metadata columns are now consistently set as factors.

## bcbioRNASeq 0.1.1 (2017-10-26)

- Added support for coloring of multiple interesting groups in quality control
  plots.

## bcbioRNASeq 0.1.0 (2017-10-23)

- Updated version and author information to match the F1000 Research
  workflow.
- Added an `f1000v1` branch containing the reproducible code used to generate
  the figures in our workflow.
- Modified `plotMA` to support vertical or horizontal layout return. Also
  added an argument to remove the color legend, which is typically not that
  informative.
- Added custom color palette support to the quality control functions.
- Upgrading from `bcbioRNADataSet` (< 0.1.0) to `bcbioRNASeq` class object is
  now possible using `as` coercion method.
- Object oriented methods are now restricted to use `bcbioRNASeq` object.
  Legacy `bcbioRNADataSet` objects must be upgraded to `bcbioRNASeq` class.

## bcbioRNASeq 0.0.28 (2017-10-17)

- Added support for output of unstructured data inside `bcbioRNASeq` S4 object
  using `flatFiles` function.
- Added `bcbioRNASeq` method support for `annotable` generic.

## bcbioRNASeq 0.0.27 (2017-10-10)

- Renamed `bcbioRNADataSet` S4 class to `bcbioRNASeq`. This matches the
  naming conventions in the [bcbioSingleCell][] package.
- Renamed `loadRNASeqRun` to simply `loadRNASeq`.
- Switched `loadRNASeq` from using S4 dispatch to a standard function.
- Added a parameter argument to `loadRNASeq` that enables request of a
  specific Ensembl release version for gene annotations.
- Renamed `interestingGroup` argument in quality control functions to
  `interestingGroups` for better consistency.
- Improved handling of sample metrics in `plotPCACovariates`.
- Added functional analysis R Markdown template.
- Offloaded some core functionality shared between [bcbioRNASeq][] and
  [bcbioSingleCell][] to the [basejump][] package. This included some code to
  handle sample metadata YAML and file loading. This helps provide a consistent
  experience across both packages.

## bcbioRNASeq 0.0.26 (2017-09-09)

- Renamed package from bcbioRnaseq to bcbioRNASeq.
- Improved website appearance.
- Added [viridis][] color palette support to quality control functions.
- Improved subset operations on `bcbioRNADataSet` object.
- Fixed setup chunk loading of `bcbioRNADataSet` in differential expression
  [R Markdown][] template.

## bcbioRNASeq 0.0.25 (2017-08-11)

- Added S4 methods support for plots, allowing the user to use either
  `bcbioRNADataSet` or a metrics `data.frame` and manual `interesting_group`
  declaration for visualization.
- Migrated function and variable names from `snake_case` to `camelCase`.
- Offloaded small RNA functionality to a separate package named
  [bcbioSmallRNA][].

## bcbioRNASeq 0.0.24 (2017-07-13)

- Reworked [R Markdown][] templates to improve YAML defaults and add more
  comments.
- Modified default path variables in `setup.R` to use `*_dir` instead of
  `*_out`.
- Updated NEWS file to use [Markdown][] syntax.

## bcbioRNASeq 0.0.23 (2017-07-03)

- Slotted `DESeqDataSet` using `design = formula(~1)` for quality control. This
  enables automatic generation of `rlog` and `vst` transformed counts.
- Documentation fixes and website updates.
- Renamed S4 class from `bcbioRnaDataSet` to `bcbioRNADataSet` (case sensitive).
- Adjusted the number of exported functions.

## bcbioRNASeq 0.0.22 (2017-06-21)

- Added [testthat][] checking with [lintr][].
- Initial setup of code coverage using [covr][].

## bcbioRNASeq 0.0.21 (2017-06-16)

- Prepared draft of [F1000][] workflow document.

## bcbioRNASeq 0.0.20 (2017-06-09)

- Added Travis CI support for automatic rendering of quality control report.

## bcbioRNASeq 0.0.19 (2017-06-07)

- `bcbioRnaDataSet` S4 definition updates.
- Updates to `plot_pca` and gene-level heatmaps.

## bcbioRNASeq 0.0.18 (2017-05-24)

- Simplified count pooling functions.

## bcbioRNASeq 0.0.17 (2017-05-21)

- Reduced number of exports and improved documentation.

## bcbioRNASeq 0.0.16 (2017-05-18)

- Draft migration of [bcbio][] run object into S4 `bcbioRnaDataSet`.
- Created a new variant of `load_run` that saves to S4 object instead of list.

## bcbioRNASeq 0.0.15 (2017-05-15)

- Reworked and re-organized internal functions.

## bcbioRNASeq 0.0.14 (2017-05-10)

- Defaulted to loading run using project summary YAML file.
- Initial commit of [R Markdown][] templates (e.g. quality control).
- Added support for dynamic file downloads from [HBC][] website.
- Draft build of website using `pkgdown::build_site`.

## bcbioRNASeq 0.0.13 (2017-05-08)

- Improved [RDAVIDWebService][] utility functions to work with [dplyr][] 0.6.0.

## bcbioRNASeq 0.0.12 (2017-05-01)

- Reworked metadata and summary metrics functions to obtain information from
  `project-summary.yaml` saved in the final run directory.

## bcbioRNASeq 0.0.11 (2017-04-27)

- Reduced number of depdencies.
- Initial commit of modified volcano plot from [CHBUtils][] package.
- Internal code updates for upcoming [dplyr][] 0.6.0/[tidyeval][] update.
- Updated [Ensembl][] [biomaRt][] annotations to use live site, currently
  release 88.

## bcbioRNASeq 0.0.10 (2017-04-19)

- Renamed `import_*` functions to `read_*`.

## bcbioRNASeq 0.0.9 (2017-04-13)

- Consolidated NAMESPACE imports.
- Defaulted to writing count matrices with gzip compression, to save disk space.

## bcbioRNASeq 0.0.8 (2017-04-12)

- Renamed internal parameters for better readability.
- Improved documentation and consolidate functions by group.

## bcbioRNASeq 0.0.7 (2017-04-10)

- NAMESPACE simplification using [basejump][] package.

## bcbioRNASeq 0.0.6 (2017-04-07)

- Reworked handling of plots and tables during knits.

## bcbioRNASeq 0.0.5 (2017-04-06)

- Initial commit of differential expression and gene set enrichment functions.

## bcbioRNASeq 0.0.4 (2017-04-04)

- Added [bcbio][] object integrity checks.
- Improved detection and handling of lane split samples.

## bcbioRNASeq 0.0.3 (2017-03-31)

- Reworked functions to utilize [bcbio][] list object.

## bcbioRNASeq 0.0.2 (2017-03-28)

- Added plotting functions.

## bcbioRNASeq 0.0.1 (2017-03-22)

- Start of package development.
- Initial draft release supporting automatic loading of [bcbio][] run data.

[acid genomics]: https://acidgenomics.com/
[annotationhub]: https://bioconductor.org/packages/AnnotationHub/
[appveyor ci]: https://ci.appveyor.com/
[basejump]: https://r.acidgenomics.com/packages/basejump/
[bcbio]: https://github.com/bcbio/bcbio-nextgen/
[bcbiobase]: https://r.acidgenomics.com/packages/bcbiobase/
[bcbiornaseq]: https://r.acidgenomics.com/packages/bcbiornaseq/
[bcbiosinglecell]: https://r.acidgenomics.com/packages/bcbiosinglecell/
[bcbiosmallrna]: https://github.com/lpantano/bcbioSmallRna/
[bioconda]: https://bioconda.github.io/
[bioconductor]: https://bioconductor.org/
[biomart]: https://bioconductor.org/packages/biomaRt/
[chbutils]: https://github.com/hbc/CHBUtils/
[clusterprofiler]: https://bioconductor.org/packages/clusterProfiler/
[covr]: https://github.com/r-lib/covr/
[degreport]: https://bioconductor.org/packages/DEGreport/
[deseq2]: https://bioconductor.org/packages/DESeq2/
[deseqanalysis]: https://r.acidgenomics.com/packages/deseqanalysis/
[dplyr]: https://dplyr.tidyverse.org/
[dropbox]: https://www.dropbox.com/
[edger]: https://bioconductor.org/packages/edgeR/
[ensdb.hsapiens.v75]: https://bioconductor.org/packages/EnsDb.Hsapiens.v75/
[ensembl]: https://www.ensembl.org/
[ensembldb]: https://bioconductor.org/packages/ensembldb/
[f1000]: https://f1000.com/
[ggplot2]: https://ggplot2.tidyverse.org/
[hbc]: https://bioinformatics.sph.harvard.edu/
[lintr]: https://github.com/r-lib/lintr/
[markdown]: https://daringfireball.net/projects/markdown/
[r markdown]: https://rmarkdown.rstudio.com/
[r]: https://www.r-project.org/
[rdavidwebservice]: https://bioconductor.org/packages/RDAVIDWebService/
[testthat]: https://github.com/r-lib/testthat/
[tidyverse]: https://www.tidyverse.org/
[tximport]: https://bioconductor.org/packages/tximport/
[viridis]: https://cran.r-project.org/package=viridis
