#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param fill Desired ggplot fill scale. Defaults to
#'   [scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [scale_fill_manual()].
#' @param flip Flip x and y axes.
#' @param headerLevel R Markdown header level.
#' @param i An integer or numeric scalar.
#' @param interestingGroups Character vector denoting groups of interest that
#'   define the samples. If left unset, defaults to `sampleName`.
#' @param object Object.
#' @param passLimit Threshold to plot pass color marker.
#' @param return Object class to return. Uses [match.arg()] internally and picks
#'   the first item in the vector by default.
#' @param title Title.
#' @param value Value to assign.
#' @param warnLimit Threshold to plot warning color marker.
#' @param withDimnames A `logical`, indicating whether dimnames should be
#'   applied to extracted assay elements.
#' @param x Object.
#' @param ... Additional arguments.
#'
#' @return No value.
NULL
