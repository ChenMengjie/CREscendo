
annotate_peaks_with_CREs <- function(peaks.gr, CRE_gr) {

  overlaps <- findOverlaps(peaks.gr, CRE_gr)

  # Extract indices for both query and subject hits
  query_indices <- queryHits(overlaps)
  subject_indices <- subjectHits(overlaps)

  # unique_query_indices <- unique(query_indices)
  kkk_tab <- table(query_indices)
  multi_peaks <- as.numeric(names(kkk_tab)[which(kkk_tab>=2)])
  mtone_peaks <- as.numeric(names(kkk_tab)[which(kkk_tab>=1)])

  return(list(overlaps = overlaps, multi_CREs_peaks = multi_peaks, mtone_CREs_peaks = mtone_peaks))
}



perform_chisq_on_CREs <- function(all_peaks_summary, peaks_annotated_with_CREs, CRE_gr, multiCREs){

  annotation_using_ENCODE_CRE <- vector("list", length(all_peaks_summary))

  overlap_info <- peaks_annotated_with_CREs$overlaps
  query_id <- queryHits(overlap_info)
  subject_id <- subjectHits(overlap_info)

  chisq_test_on_CREs <- list(NULL)

  for(i in 1:length(all_peaks_summary)){
    if(multiCREs){
      peakID <- peaks_annotated_with_CREs$multi_CREs_peaks[i]
    } else {
      peakID <- peaks_annotated_with_CREs$mtone_CREs_peaks[i]
    }
    CREID <- subject_id[query_id %in%peakID]
    peaksENCODE <- CRE_gr[CREID]

    # Extract start and end positions of overlapping regions
    starts <- start(peaksENCODE)
    peaksENCODE <- peaksENCODE[order(starts)]
    CREID <- CREID[order(starts)]

    starts <- start(peaksENCODE)
    ends <- end(peaksENCODE)
    if(nrow(all_peaks_summary[[i]]) > 0){
      # Initialize group_info with zeros
      group_info <- integer(nrow(all_peaks_summary[[i]]))

      # Assign each position in all_peaks_summary[[i]] to a group based on overlap
      for (j in seq_along(starts)) {
        # Vectorized operation to assign groups
        group_info[all_peaks_summary[[i]]$pos >= starts[j] & all_peaks_summary[[i]]$pos <= ends[j]] <- j
      }

      sub_data <- all_peaks_summary[[i]][, -1]

      if(nrow(sub_data) == 0 | sum(sub_data[, 1]) == 0 |sum(sub_data[, 2])  == 0 |  length(unique(group_info)) == 1){
        chisq_test_on_CREs[[i]] <- NA
      } else {

        #### test ####
        merged_data <- apply(sub_data, 2, function(x){
          tapply(x, group_info, sum)
        })

        freq_vec <- apply(merged_data, 2, function(x){x/sum(x)})
        merged_data[is.na(merged_data)] <- 0

        flag <- which(merged_data[, 1]>0)
        kkk1 <- contigency(merged_data[flag, 2], merged_data[flag, 1])
        effective_df1 <- length(kkk1)
        chisq_stat1 <- sum(kkk1^2)
        contribution1 <- 100*kkk1^2/chisq_stat1

        flag <- which(merged_data[, 2]>0)
        kkk2 <- contigency(merged_data[flag, 1], merged_data[flag, 2])
        effective_df2 <- length(kkk2)
        chisq_stat2 <- sum(kkk2^2)
        contribution2 <- 100*kkk2^2/chisq_stat2

        chisq_summary <- c(n_1 = sum(merged_data[, 1]), n_2 = sum(merged_data[, 2]),
                           chisq_stat1 = chisq_stat1, df1 = effective_df1, pval1 = 1-pchisq(chisq_stat1, df = effective_df1-1),
                           chisq_stat2 = chisq_stat2, df2 = effective_df2, pval2 = 1-pchisq(chisq_stat2, df = effective_df2-1))

        chisq_test_on_CREs[[i]] <- list(peakID = peakID, CREID = CREID,
                                        CRE_counts = merged_data,
                                        group_info = group_info,
                                        freq_vec = freq_vec,
                                        chisq_summary = chisq_summary,
                                        contribution1 = contribution1,
                                        contribution2 = contribution2)
      }
    } else {
      chisq_test_on_CREs[[i]] <- NA
    }
    names(chisq_test_on_CREs)[[i]] <- names(all_peaks_summary)[[i]]
  }

  return(chisq_test_on_CREs)
}


summarize_chisq_test_on_CREs <- function(chisq_test_on_CREs){

  chi_test_summary <- lapply(chisq_test_on_CREs, function(x){
    if(length(x)==1){
      return(c(n_1 = NA, n_2 = NA, chisq_stat1 = NA, df1 = NA, pval1 = NA,
               chisq_stat2 = NA, df2 = NA, pval2 = NA))
    } else {
      return(x$chisq_summary)
    }})

  chi_test_summary <- do.call(rbind, chi_test_summary)
  rownames(chi_test_summary) <- names(chisq_test_on_CREs)

  ####### chisq_test summary on CREs ########
  CRE_summary <- lapply(chisq_test_on_CREs, function(x){
    if(length(x)==1){
      return(c(max_dif = NA, abs_dif = NA))
    } else {
      dif <- apply(x$freq_vec, 1, function(x){
        abs(x[1]-x[2])/max(x[1], x[2])
      })
      abs_dif <- apply(x$freq_vec, 1, function(x){
        abs(x[1]-x[2])
      })

      return(c(max_dif =  max(dif), abs_dif = max(abs_dif),
             max_contribution1 = max(x$contribution1),
             max_contribution2 = max(x$contribution2)))
    }})

  CRE_summary <- do.call(rbind, CRE_summary)
  chi_test_summary <- cbind(chi_test_summary, CRE_summary)
  chi_test_summary <- as.data.frame(chi_test_summary)

  chi_test_summary <- chi_test_summary[!is.na(chi_test_summary[, 1]), ]
  chi_test_summary$fdr1 <- p.adjust(chi_test_summary$pval1, method = "fdr")
  chi_test_summary$fdr2 <- p.adjust(chi_test_summary$pval2, method = "fdr")
  fdrs <- cbind(chi_test_summary$fdr1, chi_test_summary$fdr2)
  chi_test_summary$max_refID <- apply(fdrs, 1, function(x){which.max(x)})
  chi_test_summary$chisq_stat <- ifelse(chi_test_summary$max_refID == 1,
                                        chi_test_summary$chisq_stat1, chi_test_summary$chisq_stat2)
  chi_test_summary$fdr <- ifelse(chi_test_summary$max_refID == 1,
                                        chi_test_summary$fdr1, chi_test_summary$fdr2)
  chi_test_summary$max_contribution <- ifelse(chi_test_summary$max_refID == 1,
                                 chi_test_summary$max_contribution1, chi_test_summary$max_contribution2)
  chi_test_summary <- chi_test_summary[order(chi_test_summary$chisq_stat, decreasing = T), ]
  chi_test_summary <- chi_test_summary[, c("chisq_stat", "fdr", "max_dif", "abs_dif", "max_contribution", "max_refID",
                                           "n_1", "n_2",
                                           "chisq_stat1", "df1", "pval1", "fdr1",
                                           "chisq_stat2", "df2", "pval2", "fdr2")]
  return(chi_test_summary)

}


chisq_test_on_peaks_using_CRE_annotations <- function(all_peaks_summary,
                                                      peaks_annotated_with_CREs, CRE_gr,
                                                      celltypes = NULL, multiCREs){

  if(ncol(all_peaks_summary[[1]]) == 3){
    chisq_test_on_CREs <- perform_chisq_on_CREs(all_peaks_summary, peaks_annotated_with_CREs, CRE_gr, multiCREs)
    chisq_test_on_CREs_summary <- summarize_chisq_test_on_CREs(chisq_test_on_CREs)
    res <- list(chisq_test_on_CREs = chisq_test_on_CREs, chisq_test_on_CREs_summary = chisq_test_on_CREs_summary)

  } else if (ncol(all_peaks_summary[[1]]) > 3){

    if(is.null(celltypes) | length(celltypes) > 2){
      stop("More than cell types are available. Please specify one or two cell types.")
    } else if(length(celltypes) == 2){
      if(!all(celltypes%in%colnames(all_peaks_summary[[1]])))
        stop("Please specify the correct cell types.")
      else
        message("Two cell types are specified. Comparing peaks between the two cell types.")
    } else if(length(celltypes) == 1){
      if(!celltypes%in%colnames(all_peaks_summary[[1]]))
        stop("Please specify the correct cell types.")
      else
        message("One cell type is specified. Comparing peaks between the specified cell type and the rest.")
    }

    if(length(celltypes) == 2){
      peaks_summary <- lapply(all_peaks_summary, function(x){
        x <- x[, c("pos", celltypes)]
        return(x)
      })
    }
    if(length(celltypes) == 1){
      peaks_summary <- lapply(all_peaks_summary, function(x){
        y <- x[, !colnames(x)%in%c("pos", celltypes)]
        sum_y <- rowSums(y)
        x <- cbind(x[, "pos", drop = F], x[, celltypes, drop = F], rest = sum_y)
        colnames(x) <- c("pos", celltypes, "rest")
        return(x)
      })
    }
    chisq_test_on_CREs <- perform_chisq_on_CREs(peaks_summary, peaks_annotated_with_CREs, CRE_gr, multiCREs)
    chisq_test_on_CREs_summary <- summarize_chisq_test_on_CREs(chisq_test_on_CREs)
    res <- list(chisq_test_on_CREs = chisq_test_on_CREs, chisq_test_on_CREs_summary = chisq_test_on_CREs_summary)
  }
  return(res)
}

annotate_chisq_test_results <- function(chisq_test_on_peaks, annot_peaks, CRE_gr) {
  # Extract summary and CREs test results
  chisq_test_on_CREs <- chisq_test_on_peaks$chisq_test_on_CREs
  chisq_test_on_CREs_summary <- chisq_test_on_peaks$chisq_test_on_CREs_summary

  # Get row names from the summary
  row_names <- rownames(chisq_test_on_CREs_summary)

  # Use lapply to iterate over row names
  results <- lapply(row_names, function(which_peak) {
    CRE_info <- CRE_gr[chisq_test_on_CREs[[which_peak]]$CREID]
    CRE_type <- paste(CRE_info$type, collapse = ";")
    pID <- chisq_test_on_CREs[[which_peak]]$peakID

    return(list(CRE_type = CRE_type, pID = pID))
  })

  # Unlist the results for CRE_type and pID separately
  CRE_types <- sapply(results, function(x) x$CRE_type)
  pIDs <- sapply(results, function(x) x$pID)

  # Combine the results into a new data frame
  chisq_test_summary <- data.frame(
    chisq_test_on_CREs_summary[, c("chisq_stat", "fdr", "max_contribution", "max_refID", "max_dif", "abs_dif", "n_1", "n_2", "df1")],
    CRE_type = CRE_types,
    annot_peaks[as.character(pIDs), ]
  )

  return(list(chisq_test_on_CREs = chisq_test_on_CREs, chisq_test_on_CREs_summary = chisq_test_summary))
}



annotated_tests_on_peaks <- function(all_peaks_summary, annot_peaks, peaks_annotated_with_CREs, CRE_gr, celltypes = NULL, multiCREs){

  chisq_test_on_peaks <- chisq_test_on_peaks_using_CRE_annotations(all_peaks_summary,
                                                      peaks_annotated_with_CREs, CRE_gr,
                                                      celltypes, multiCREs)
  chisq_test_on_peaks <- annotate_chisq_test_results(chisq_test_on_peaks, annot_peaks, CRE_gr)

  return(chisq_test_on_peaks)

}


extract_peak_CREs_cleavage <- function(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary, CRE_gr, which_peak){

  CREID <- chisq_test_on_CREs[[which_peak]]$CREID
  CRE_info <- CRE_gr[CREID]
  freq_info <- chisq_test_on_CREs[[which_peak]]$freq_vec
  cleavage <- all_peaks_summary[[which_peak]]
  annotation <- chisq_test_summary[which_peak, ]
  freq <- chisq_test_on_CREs[[which_peak]]$CRE_counts[-1, ]
  freq <- apply(freq, 2, function(x){
    if(sum(x) == 0){
      rep(0, length(x))
    }
    else {
      x/sum(x)
    }
    })
  merged_fr <- data.frame(CRE_id = paste0("CRE", rownames(freq)), freq)
  colnames(merged_fr) <- c("CRE_id", colnames(freq_info))
  merge_summary <- melt(merged_fr, id = "CRE_id")
  colnames(merge_summary)[2:3] <- c("celltype", "frequency")

  return(list(CRE_info = CRE_info, freq_info = freq_info,
              freq_CRE = merge_summary, cleavage = cleavage, annotation = annotation))
}



extract_gene_CREs_cleavage <- function(chisq_test_on_CREs, chisq_test_summary, all_peaks_summary, CRE_gr, which_gene){

  sub_annot <- chisq_test_summary[chisq_test_summary$SYMBOL%in%which_gene, ]
  which_peaks <- rownames(sub_annot)[order(rownames(sub_annot))]

  cleavages <- all_peaks_summary[which_peaks]
  cleavages <- do.call(rbind, cleavages)

  CRE_count_info <- lapply(seq_along(chisq_test_on_CREs[which_peaks]), function(index) {
    x <- chisq_test_on_CREs[[which_peaks[index]]]
    if("0"%in%rownames(x$CRE_counts)){
      freq <- x$CRE_counts[-1, ]
    } else {
      freq <- x$CRE_counts
    }
    if(is.null(nrow(freq))){
      freq_CRE_each_peak <- data.frame(peakID = index, matrix(1, ncol = length(freq)))  # Adding peakID column
      colnames(freq_CRE_each_peak)[-1] <- names(freq)
    } else {

      freq <- apply(freq, 2, function(x){
       if(sum(x) == 0){
         rep(0, length(x))
        }
        else {
         x/sum(x)
        }
     })
      freq_CRE_each_peak <- data.frame(peakID = index, freq)  # Adding peakID column
      colnames(freq_CRE_each_peak)[-1] <- colnames(freq)
    }

    return(freq_CRE_each_peak)
  })

  CRE_count_info <- do.call(rbind, CRE_count_info)

  merged_fr <- data.frame(CRE_id = paste0("CRE", 1:nrow(CRE_count_info)), CRE_count_info)
  colnames(merged_fr) <- c("CRE_id", colnames(CRE_count_info))
  merge_summary <- melt(merged_fr, id = c("CRE_id", "peakID"))
  colnames(merge_summary)[3:4] <- c("celltype", "frequency")

  CRE_ID <- lapply(chisq_test_on_CREs[which_peaks], function(x) x$CREID)
  CRE_ID <- unlist(CRE_ID)
  CRE_info <- CRE_gr[CRE_ID]



  return(list(CRE_info = CRE_info, CRE_freq_info = merge_summary,
              cleavages = cleavages, peak_annotation = sub_annot[which_peaks, ]))

}


