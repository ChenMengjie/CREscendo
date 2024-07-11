
#' Create a CREscendo Object
#'
#' This constructor function initializes a CREscendo object containing various components
#' related to cis-regulatory element (CRE) analysis. The object consolidates annotated peaks,
#' counts of cleavage sites, CRE annotations, and peaks annotated with CREs.
#'
#' @param annotated_peaks A data frame containing information about the annotated peaks.
#' @param cleavage_sites_counts A list of matrix with counts of cleavage sites across different peak regions.
#' @param CRE_annotations Data frame detailing annotations of CREs, including type, location, etc.
#' @param peaks_annotated_with_CREs A GenomicRange that links peaks to specific CRE annotations.
#' @param multiCREs A logical value indicating whether to include peaks with multiple CREs (TRUE) or only peaks with one CRE (FALSE).
#' @return An object of class 'CREscendo' containing the specified datasets, structured for further analysis.
#'
#' @export
CreateCREscendoObject <- function(annotated_peaks, cleavage_sites_counts,
                                  CRE_annotations, peaks_annotated_with_CREs, multiCREs) {
  structure(list(
    annotated_peaks = annotated_peaks,
    cleavage_sites_counts = cleavage_sites_counts,
    CRE_annotations = CRE_annotations,
    peaks_annotated_with_CREs = peaks_annotated_with_CREs,
    multiCREs = multiCREs
  ), class = "CREscendo")
}

#' Prepare a CREscendo Object for Analysis
#'
#' This function prepares and constructs a CREscendo object by processing genomic peaks,
#' associating them with cis-regulatory elements (CREs), annotating these peaks,
#' and compiling counts of cleavage sites for specified cell types. It integrates multiple
#' genomic and annotation resources to create a comprehensive data object for further analysis.
#'
#' @param peaks_gr GenomicRanges object containing genomic peak locations.
#' @param tabix_file Path to a Tabix-indexed file containing fragments relevant to the peaks.
#' @param cell_annotation Data frame containing cell annotations with columns for cellname and celltype.
#' @param celltype_list Optional; a vector of cell types to include in the analysis. If NULL,
#'        all unique cell types from cell_annotation will be used.
#' @param CRE_gr GenomicRanges object for cis-regulatory elements.
#' @param tssRegion Numeric vector of length 2 specifying the region around transcription start sites
#'        to consider for annotation (defaults to c(-1000, 1000)).
#' @param TxDb A TxDb object containing transcript annotations.
#' @param annoDb Character string, the Bioconductor annotation package for organism specific data, default
#'        is "org.Hs.eg.db" for human.
#' @param multiCREs Logical value indicating whether to include peaks with multiple CREs (TRUE) or all peaks (FALSE).
#'
#' @return A CREscendo object containing annotated peaks, CRE associations, and cleavage sites counts
#'         tailored for further analysis.
#' @importFrom ChIPseeker annotatePeak
#' @export
PrepareCREscendoObject <- function(peaks_gr, tabix_file, cell_annotation, celltype_list,
                                   CRE_gr, tssRegion=c(-1000, 1000), TxDb, annoDb="org.Hs.eg.db", multiCREs = TRUE){

  peaks_annotated_with_CREs <- annotate_peaks_with_CREs(peaks_gr, CRE_gr)
  if(multiCREs){
    seleced_peakIDs <- peaks_annotated_with_CREs$multi_CREs_peaks
    selected_peaks <- peaks_gr[seleced_peakIDs]
  } else {
    seleced_peakIDs <- peaks_annotated_with_CREs$mtone_CREs_peaks
    selected_peaks <- peaks_gr[seleced_peakIDs]
  }
  bed_annot <- ChIPseeker::annotatePeak(selected_peaks, tssRegion=tssRegion, TxDb=TxDb, annoDb=annoDb)
  annot_peaks <- as.data.frame(bed_annot)
  cellnames <- cell_annotation$cellname
  celltypes <- cell_annotation$celltype
  if(is.null(celltype_list)){
    celltype_list <- unique(celltypes)
  }
  cleavage_sites_counts <- extracted_cleavage_for_a_list_of_peaks_for_celltypes(annot_peaks, tabix_file, cellnames, celltypes, celltype_list)
  rownames(annot_peaks) <- seleced_peakIDs
  res <- CreateCREscendoObject(annot_peaks, cleavage_sites_counts, CRE_gr, peaks_annotated_with_CREs, multiCREs)
  return(res)

}

#' Create a CREscendoTested Object
#'
#' This constructor function initializes a CREscendoTested object, which encapsulates the results
#' and relevant data from CRE association tests. This object includes information about cell types,
#' chi-squared test results, annotated peaks, cleavage sites counts, CRE annotations, and
#' peaks associated with CREs.
#'
#' @param CREscendo A CREscendo object containing initial data such as annotated peaks and
#' cleavage sites counts.
#' @param celltypes A vector of cell types for which tests were conducted. This parameter also
#' automatically includes a negation of each cell type (e.g., "non-Type1" if "Type1" is specified)
#' if only one cell type is provided.
#' @param annotated_tests_on_peaks A list containing results of tests on peaks, including chi-squared
#' test results and summaries.
#'
#' @return A CREscendoTested object that contains detailed results and annotations from CRE testing
#' processes, structured for further analysis or visualization.
#'
#' @export
CreateCREscendoTestedObject <- function(CREscendo, celltypes, annotated_tests_on_peaks) {
  if(length(celltypes) == 1){
    celltypes <- c(celltypes, paste0("non-", celltypes))
  }
  structure(list(
    celltypes = celltypes,
    chisq_test_on_CREs = annotated_tests_on_peaks$chisq_test_on_CREs,
    chisq_test_summary = annotated_tests_on_peaks$chisq_test_on_CREs_summary,
    annotated_peaks = CREscendo$annotated_peaks,
    cleavage_sites_counts = CREscendo$cleavage_sites_counts,
    CRE_annotations = CREscendo$CRE_annotations,
    peaks_annotated_with_CREs = CREscendo$peaks_annotated_with_CREs
  ), class = "CREscendoTested")
}

#' Print method for CREscendo objects
#'
#' This function provides a concise summary of a CREscendo object when it is printed to the console.
#' It displays the number of peaks that overlap with more than one cis-regulatory element (CRE)
#' and lists the cell types for which fragment cleavage sites have been extracted.
#'
#' @param x A CREscendo object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The output includes:
#' - The number of peaks that have overlaps with multiple CREs.
#' - A list of cell types for which fragment cleavage site data is available.
#' @export
#' @method print CREscendo
print.CREscendo <- function(x, ...) {
  if(!inherits(x, "CREscendo"))
    stop("The input object must be a CREscendo object.")
  cat("CREscendo Object:\n")
  cat(nrow(x$annotated_peaks), " peaks overlapped with more than 1 cis-regularory element (CRE).\n")
  cat("Fragment cleavage sites extracted for the follow cell types: \n")
  print(colnames(x$cleavage_sites_counts[[1]])[-1])
}


#' Test CRE Association with specified cell types within a CREscendo Object
#'
#' This function performs chi-squared tests on cis-regulatory elements (CREs) associated with peaks
#' within a CREscendo object. It evaluates differential CRE activities between two cell types.
#' This is intended to facilitate the analysis of CRE dynamics and their potential regulatory impacts.
#'
#' @param x A CREscendo object containing the necessary datasets for CRE analysis.
#' @param celltypes One or two cell types. Must be cell types included in the CREscendo object.
#' If one cell type is specified, the comparison will be made between the specified cell type and all other cell types.
#' If two cell types are specified, the comparison will be made between the two specified cell types.
#' @param ... Additional arguments passed to lower-level functions, such as plotting or statistical functions.
#' @return Returns a CREscendoTestedObject containing the results of the CRE tests,
#' which include statistical summaries and any relevant annotations for further analysis.
#' @export
#' @method CREtest CREscendo
CREtest.CREscendo <- function(x, ..., celltypes = NULL) {

  if(!inherits(x, "CREscendo"))
    stop("The input object must be a CREscendo object.")
  # Extract the data from the CREscendo object
  annot_peaks <- x$annotated_peaks
  all_peaks_summary <- x$cleavage_sites_counts
  CRE_gr <- x$CRE_annotations
  peaks_annotated_with_CREs <- x$peaks_annotated_with_CREs
  multiCREs <- x$multiCREs
  res <- annotated_tests_on_peaks(all_peaks_summary, annot_peaks, peaks_annotated_with_CREs, CRE_gr, celltypes = celltypes, multiCREs = multiCREs)
  y <- CreateCREscendoTestedObject(x, celltypes, res)
  return(y)

}

#' Summarize the Results of CREscendo Tested Object
#'
#' This function retrieves the summary of chi-squared test results from a `CREscendoTested` object.
#' It is specifically used to extract and return the chi-squared test summary component,
#' allowing for quick access to the statistical analysis results.
#'
#' @param x A `CREscendoTested` object containing results from CRE association tests.
#' @param ... Additional arguments passed to lower-level functions, such as plotting or statistical functions.
#' @return Returns a data frame containing
#' the summary of chi-squared tests conducted on CREs.
#' @export
#' @method Summarize CREscendoTested
Summarize.CREscendoTested <- function(x, ...){
  if(!inherits(x, "CREscendoTested"))
    stop("The input object must be a CREscendoTested object.")
  return(x$chisq_test_summary)
}


#' Filter the Results of CREscendo Tested Object
#'
#' This function retrieves the summary of chi-squared test results from a `CREscendoTested` object.
#' It is specifically used to extract and return high confident results from chi-squared test summary component,
#' allowing for prioritized downstream analysis.
#'
#' @param x A `CREscendoTested` object containing results from CRE association tests.
#' @param fdr_cutoff A numeric value specifying the fdr adjusted p-value threshold for filtering the results. The default is 0.05.
#' @param coverage_cutoff A numeric value specifying the cut-off of the number of fragments across all cells within the same type for each peak region for filtering the results. The default is 100.
#' @param abs_diff_cutoff A numeric value specifying the cut-off of the absolute difference in the CRE frequencies between two cell types. The default is 0.2.
#' @param ... Additional arguments passed to lower-level functions, such as plotting or statistical functions.
#' @return Returns a data frame containing
#' the summary of chi-squared tests conducted on CREs after applying several filters.
#' The test results are sorted by their chi-squared test statistics in descending order.
#' @export
#' @method Filter CREscendoTested
Filter.CREscendoTested <- function(x, ...,
                                   fdr_cutoff = 0.05, coverage_cutoff = 100, abs_diff_cutoff = 0.2){
  if(!inherits(x, "CREscendoTested"))
    stop("The input object must be a CREscendoTested object.")
  aa <- x$chisq_test_summary
  aa <- aa[aa$fdr <= fdr_cutoff & aa$n_1 >= coverage_cutoff & aa$n_2 >= coverage_cutoff & aa$abs_dif >= abs_diff_cutoff, ]
  return(aa)
}


#' Visualize cleavage sites of different cell types for a specific peak in a CREscendoTested Object
#'
#' This function provides a visualization of a specific peak within a CREscendoTested object. It extracts relevant data for the specified peak, including cleavage sites,
#' CRE annotation and CRE frequencies within tested cell types.
#'
#' @param x A `CREscendoTested` object containing results from differential CRE tests, including
#' statistical summaries and annotations.
#' @param which_peak A character string or numeric identifier specifying the peak for which
#' visualization is desired. This should correspond to a peak ID contained within the `CREscendoTested` object.
#' @param which_gene A character string specifying the gene symbol for which
#' visualization is desired. This should correspond to a gene symbol contained within the `CREscendoTested` object.
#' @param ... Additional arguments passed to lower-level functions, such as plotting or statistical functions.
#' @return A plot object representing the results for the specified peak. This allows
#' for further manipulation or saving of the plot outside the function.
#'
#' @details The function first checks that the input is a `CREscendoTested` object. If not,
#' it stops with an error message. It then extracts the data necessary for visualization from
#' the object, including chi-squared test results and cleavage sites counts for the specified peak or the peaks within a specific gene.
#' Finally, it calls internal functions to generate the visualization, which is then returned.
#'
#' @export
#' @method Visualize CREscendoTested
Visualize.CREscendoTested <- function(x, ..., which_peak = NULL, which_gene = NULL){
  if(!inherits(x, "CREscendoTested"))
    stop("The input object must be a CREscendoTested object.")

  if(is.null(which_peak) & is.null(which_gene))
    stop("Please specify a peak ID or a gene name for visualization.")

  if(!is.null(which_peak) & !is.null(which_gene))
    stop("Please specify a peak ID or a gene name for visualization.")

  if(!is.null(which_peak) & is.null(which_gene) ){

    if(!which_peak%in%rownames(x$chisq_test_summary))
      stop("The specified peak ID is not found in the CREscendoTested object.")
    else {
      chisq_test_on_CREs <- x$chisq_test_on_CREs
      chisq_test_summary <- x$chisq_test_summary
      all_peaks_summary <- x$cleavage_sites_counts
      CRE_gr <- x$CRE_annotations
      toplot <- extract_peak_CREs_cleavage(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary,
                                        CRE_gr, which_peak)
      pp <- display_a_peak(toplot, celltypes = x$celltypes)
    }
  }

  if(is.null(which_peak) & !is.null(which_gene)){

    if(!which_gene%in%x$chisq_test_summary$SYMBOL)
      stop("The specified gene name is not found in the CREscendoTested object.")
    else {
      chisq_test_on_CREs <- x$chisq_test_on_CREs
      chisq_test_summary <- x$chisq_test_summary
      all_peaks_summary <- x$cleavage_sites_counts
      CRE_gr <- x$CRE_annotations
      gene_toplot <- extract_gene_CREs_cleavage(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary,
                                        CRE_gr, which_gene)
      pp <- display_a_gene(gene_toplot, celltypes = x$celltypes)
    }
  }

  return(pp)
}




#' Extract data of a specific peak or a specific gene from a CREscendoTested object
#'
#' This function extracts relevant data from a CREscendoTested object, including cleavage sites, CRE annotation, and CRE frequencies within tested cell types.
#'
#' @param x A CREscendoTested object containing results from differential CRE tests.
#' @param which_peak A character string or numeric identifier specifying the peak for visualization.
#' @param which_gene A character string specifying the gene symbol for which visualization is desired.
#' @param plot_contribution A logical value indicating whether to plot the contribution of each CRE to the Chi-square statistics.
#' @param ... Additional arguments passed to lower-level functions, such as plotting or statistical functions.
#' @return A list containing extracted data for further analysis.
#' @export
#' @details The function first checks that the input is a `CREscendoTested` object. If not,
#' it stops with an error message. It then extracts the data for the specified peak or the peaks within a specific gene from
#' the object, including chi-squared test results and cleavage sites counts. It then calls internal functions to quantify contribution of each CRE to the Chi-square statistics.
#'
#' @method Extract CREscendoTested
Extract.CREscendoTested <- function(x, ..., which_peak = NULL, which_gene = NULL, plot_contribution = TRUE){
  if(!inherits(x, "CREscendoTested"))
    stop("The input object must be a CREscendoTested object.")

  if(is.null(which_peak) & is.null(which_gene))
    stop("Please specify a peak ID or a gene name for extraction.")

  if(!is.null(which_peak) & !is.null(which_gene))
    stop("Please specify a peak ID or a gene name for extraction.")

  if(!is.null(which_peak) & is.null(which_gene) ){

    if(!which_peak%in%rownames(x$chisq_test_summary))
      stop("The specified peak ID is not found in the CREscendoTested object.")
    else {
      chisq_test_on_CREs <- x$chisq_test_on_CREs
      chisq_test_summary <- x$chisq_test_summary
      all_peaks_summary <- x$cleavage_sites_counts
      CRE_gr <- x$CRE_annotations
      extracted <- extract_peak_CREs_cleavage(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary,
                                          CRE_gr, which_peak)
      refID <- unlist(extracted$annotation[1, "max_refID"])
      if(refID == 1){
        extracted$CREcontribution <- chisq_test_on_CREs[[which_peak]]$contribution1
      } else {
        extracted$CREcontribution <- chisq_test_on_CREs[[which_peak]]$contribution2
      }
      names(extracted$CREcontribution) <- paste0("CRE", names(extracted$CREcontribution))
      names(extracted$CREcontribution)[1] <- "nonCRE"
    }
    if(plot_contribution){

      data <- data.frame(x = 1:length(extracted$CREcontribution), y = 1, value = extracted$CREcontribution)  # Replace with your own data
      # Create a scatter plot with circle size representing values
      kk <- ggplot(data, aes(x = x, y = y, size = value)) +
        geom_point(shape = 21, fill = "darkred", color = "black") +
        geom_text(aes(label = rownames(data)), size = 3, color = "black", nudge_y = 0.5) +
        scale_size_continuous(range = c(1, 15), name = "Contribution") +  # Adjust the range of circle sizes
        labs(title = which_peak, x = "CRE", y = "") + ylim(0, 2) +
        theme_void()
    }
    extracted$plot_contribution <- kk
  }
  if(is.null(which_peak) & !is.null(which_gene)){

    if(!which_gene%in%x$chisq_test_summary$SYMBOL)
      stop("The specified gene name is not found in the CREscendoTested object.")
    else {
      chisq_test_on_CREs <- x$chisq_test_on_CREs
      chisq_test_summary <- x$chisq_test_summary
      all_peaks_summary <- x$cleavage_sites_counts
      CRE_gr <- x$CRE_annotations
      extracted <- extract_gene_CREs_cleavage(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary,
                                               CRE_gr, which_gene)
    }
  }

  return(extracted)
}

