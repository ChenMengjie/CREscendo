
#### TabixOutputToDataFrame obtained from Seurat package
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  original_names = names(reads)
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0 | is.null(x = original_names)) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}


######## Extract cleavage sites for a list of peaks for celltypes
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
extracted_cleavage_for_a_list_of_peaks_for_celltypes <- function(sub_peaks, tabix.file,
                                                                cellnames, celltypes,
                                                                cell_type_list){

  # Initialize list for results
  peak_1bp_summary <- vector("list", nrow(sub_peaks))

  # Create a progress bar (update less frequently to reduce overhead)
  pb <- txtProgressBar(min = 0, max = nrow(sub_peaks), style = 3)
  interval_update <- max(1, floor(nrow(sub_peaks) / 20))  # Update progress every 5% of the total rows

  for(kkk in 1:nrow(sub_peaks)){

    aa <- sub_peaks[kkk, ]
    check_these <- GRanges(seqnames = aa$seqnames, ranges = IRanges(start = aa$start, end = aa$end))
    suppressWarnings({
      query_result <- scanTabix(file = tabix.file, param = check_these)
    })
    fragments <- TabixOutputToDataFrame(reads = query_result)
    possible_cleavage_sites <- sort(unique(c(fragments$start, fragments$end)))
    possible_cleavage_sites <- possible_cleavage_sites[possible_cleavage_sites>=aa$start &
                                                       possible_cleavage_sites<=aa$end]

    selected_intervals <- fragments[fragments$cell %in% cellnames, ]
    all_types_summary <- data.frame(pos = possible_cleavage_sites)

    results <- lapply(cell_type_list, function(celltypeID) {
      which_cell <- cellnames[celltypes %in% celltypeID]
      selected_intervals_for_type <- selected_intervals[selected_intervals$cell %in% which_cell, ]

      if(nrow(selected_intervals_for_type) == 0) {
        return(rep(NA, length(possible_cleavage_sites)))
      } else {
        cuts_summary <- table(c(selected_intervals_for_type$start, selected_intervals_for_type$end))
        return(ifelse(possible_cleavage_sites %in% names(cuts_summary), cuts_summary[as.character(possible_cleavage_sites)], 0))
      }
    })

    # Combine all cut results
    all_types_summary <- cbind(all_types_summary, do.call(cbind, results))
    all_types_summary[is.na(all_types_summary)] <- 0
    all_types_summary <- all_types_summary[rowSums(all_types_summary[, -1]) > 0, ]
    colnames(all_types_summary)[-1] <- cell_type_list
    peak_1bp_summary[[kkk]] <- all_types_summary
    names(peak_1bp_summary)[[kkk]] <- paste(aa$seqnames, aa$start, aa$end, sep = "-")
    # Update the progress bar less frequently
    if (kkk %% interval_update == 0 | kkk == nrow(sub_peaks)) {
      setTxtProgressBar(pb, kkk)
    }
  }

  close(pb)
  return(peak_1bp_summary)
}

