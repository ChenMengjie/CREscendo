% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\name{Visualize.CREscendoTested}
\alias{Visualize.CREscendoTested}
\title{Visualize cleavage sites of different cell types for a specific peak in a CREscendoTested Object}
\usage{
\method{Visualize}{CREscendoTested}(x, ..., which_peak = NULL, which_gene = NULL)
}
\arguments{
\item{x}{A `CREscendoTested` object containing results from differential CRE tests, including
statistical summaries and annotations.}

\item{...}{Additional arguments passed to lower-level functions, such as plotting or statistical functions.}

\item{which_peak}{A character string or numeric identifier specifying the peak for which
visualization is desired. This should correspond to a peak ID contained within the `CREscendoTested` object.}

\item{which_gene}{A character string specifying the gene symbol for which
visualization is desired. This should correspond to a gene symbol contained within the `CREscendoTested` object.}
}
\value{
A plot object representing the results for the specified peak. This allows
for further manipulation or saving of the plot outside the function.
}
\description{
This function provides a visualization of a specific peak within a CREscendoTested object. It extracts relevant data for the specified peak, including cleavage sites,
CRE annotation and CRE frequencies within tested cell types.
}
\details{
The function first checks that the input is a `CREscendoTested` object. If not,
it stops with an error message. It then extracts the data necessary for visualization from
the object, including chi-squared test results and cleavage sites counts for the specified peak or the peaks within a specific gene.
Finally, it calls internal functions to generate the visualization, which is then returned.
}
