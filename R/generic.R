#' Generic function for CREtest
#'
#' This function tests CRE associations in various objects.
#'
#' @param x An object to be tested.
#' @param ... Additional arguments passed to specific methods.
#'
#' @export
CREtest <- function(x, ...) {
  UseMethod("CREtest")
}


#' Generic function for Summarize
#'
#' This function provides a framework for summarizing various types of objects.
#'
#' @param x An object to summarize.
#' @param ... Additional arguments affecting the summary.
#'
#' @export
Summarize <- function(x, ...) {
  UseMethod("Summarize")
}

#' Generic function for Visualize
#'
#' This function provides a framework for visualizing various types of objects.
#'
#' @param x An object to visualize.
#' @param ... Additional arguments affecting the visualization.
#'
#' @export
Visualize <- function(x, ...) {
  UseMethod("Visualize")
}



#' Generic function for Filter
#'
#' This function provides a framework for filtering various types of objects.
#'
#' @param x An object to filter.
#' @param ... Additional arguments affecting the filtered results.
#'
#' @export
Filter <- function(x, ...) {
  UseMethod("Filter")
}





#' Generic function for Extract
#'
#' This function provides a framework for extracting information of certain peak or gene from objects.
#'
#' @param x An object to extract from.
#' @param ... Additional arguments affecting the extracting.
#'
#' @export
Extract <- function(x, ...) {
  UseMethod("Extract")
}

