# FIXME How to fix this warning?
# ❯ checking Rd cross-references ... WARNING
#   Missing link or links in documentation object 'reexports.Rd':
#     ‘[basejump]{interestingGroups<-}’ ‘[basejump]{sampleData<-}’
#
#   See section 'Cross-references' in the 'Writing R Extensions' manual.
#
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Cross_002dreferences



# 2.5 Cross-references
#
# The markup \link{foo} (usually in the combination \code{\link{foo}}) produces
# a hyperlink to the help for foo. Here foo is a topic, that is the argument of
# \alias markup in another Rd file (possibly in another package). Hyperlinks are
# supported in some of the formats to which Rd files are converted, for example
# HTML and PDF, but ignored in others, e.g. the text format.
#
# One main usage of \link is in the \seealso section of the help page, see Rd
# format.
#
# Note that whereas leading and trailing spaces are stripped when extracting a
# topic from a \alias, they are not stripped when looking up the topic of a
# \link.
#
# You can specify a link to a different topic than its name by
# \link[=dest]{name} which links to topic dest with name name. This can be used
# to refer to the documentation for S3/4 classes, for example
# \code{"\link[=abc-class]{abc}"} would be a way to refer to the documentation
# of an S4 class "abc" defined in your package, and
# \code{"\link[=terms.object]{terms}"} to the S3 "terms" class (in package
# stats). To make these easy to read in the source file,
# \code{"\linkS4class{abc}"} expands to the form given above.
#
# There are two other forms of optional argument specified as \link[pkg]{foo}
# and \link[pkg:bar]{foo} to link to the package pkg, to files foo.html and
# bar.html respectively. These are rarely needed, perhaps to refer to
# not-yet-installed packages (but there the HTML help system will resolve the
# link at run time) or in the normally undesirable event that more than one
# package offers help on a topic102 (in which case the present package has
# precedence so this is only needed to refer to other packages). They are
# currently only used in HTML help (and ignored for hyperlinks in LaTeX
# conversions of help pages), and link to the file rather than the topic (since
# there is no way to know which topics are in which files in an uninstalled
# package). The only reason to use these forms for base and recommended packages
# is to force a reference to a package that might be further down the search
# path. Because they have been frequently misused, the HTML help system looks
# for topic foo in package pkg if it does not find file foo.html.



#' @importFrom basejump export
#' @export
basejump::export

#' @importFrom basejump import
#' @export
basejump::import

#' @importFrom basejump interestingGroups
#' @export
basejump::interestingGroups

#' @rdname reexports
#' @name interestingGroups<-
#' @importFrom basejump interestingGroups<-
#' @usage NULL
#' @export
NULL

#' @importFrom basejump loadData
#' @export
basejump::loadData

#' @importFrom basejump loadRemoteData
#' @export
basejump::loadRemoteData

#' @importFrom basejump pasteURL
#' @export
basejump::pasteURL

#' @importFrom basejump sampleData
#' @export
basejump::sampleData

#' @rdname reexports
#' @name sampleData<-
#' @importFrom basejump sampleData<-
#' @usage NULL
#' @export
NULL

#' @importFrom basejump saveData
#' @export
basejump::saveData



#' @importFrom firestarter plotHeatmap
#' @export
firestarter::plotHeatmap

#' @importFrom firestarter plotQuantileHeatmap
#' @export
firestarter::plotQuantileHeatmap



#' @importFrom minimalism plotGenesDetected
#' @export
minimalism::plotGenesDetected

#' @importFrom minimalism theme_paperwhite
#' @export
minimalism::theme_paperwhite
