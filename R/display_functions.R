
make_into_tracks <- function(CRE_granges, xrange){
  xmin <- xrange[1]
  xmax <- xrange[2]
  start <- start(CRE_granges)
  end <- end(CRE_granges)

  # The palette with grey:
  cbp1 <- c( "#E69F00", "#56B4E9", "#CC79A7", "#009E73",
             "#F0E442", "#0072B2", "#D55E00",  "#999999")
  chromosome <- seqnames(CRE_granges)[1]
  to_expand <- data.frame(start = ifelse(start > xmin, start, xmin),
                          end = ifelse(end < xmax, end, xmax),
                          type = CRE_granges$type,
                          typeID = as.numeric(as.factor(CRE_granges$type)),
                          num = 1:length(CRE_granges))

  ggplot(to_expand, aes(xmin = start, xmax = end, ymin = typeID - 0.4, ymax = typeID + 0.4, fill = as.factor(typeID))) +
    geom_rect( color = "black") + xlim(xmin, xmax) +
    geom_text(data = to_expand, aes(x = (start + end) / 2, y = typeID, label = num), color = "white", size = 3)+
    #scale_y_continuous(breaks = 1:max(to_expand$typeID), labels = unique(CRE_granges$type)) +
    scale_fill_manual(values = cbp1, name = "CREs", labels = unique(CRE_granges$type)) +
    labs(x = paste0(chromosome, " position (bp)"), y = "CREs") +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),    # Remove y-axis ticks
      axis.text.y = element_blank(),     # Remove y-axis text
      axis.line.y = element_blank(),      # Remove y-axis line
      axis.ticks.x = element_blank(),    # Remove x-axis ticks
      axis.text.x = element_blank(),     # Remove x-axis text
      axis.line.x = element_blank(),      # Remove x-axis line
      axis.title.x = element_blank(),
    )
}



make_into_narrow_tracks <- function(CRE_granges, xrange){
  xmin <- xrange[1]
  xmax <- xrange[2]
  start <- start(CRE_granges)
  end <- end(CRE_granges)
  # The palette with grey:
  cbp1 <- c( "#E69F00", "#56B4E9", "#CC79A7", "#009E73",
             "#F0E442", "#0072B2", "#D55E00",  "#999999")
  chromosome <- seqnames(CRE_granges)[1]
  to_expand <- data.frame(start = ifelse(start > xmin, start, xmin),
                          end = ifelse(end < xmax, end, xmax),
                          type = CRE_granges$type,
                          typeID = as.numeric(as.factor(CRE_granges$type)),
                          num = 1:length(CRE_granges))

  ggplot(to_expand, aes(xmin = start, xmax = end, ymin = typeID -0.4, ymax = typeID + 0.4, fill = as.factor(typeID))) +
    geom_rect( color = "black") + xlim(xmin, xmax) +
   # geom_text(data = to_expand, aes(x = (start + end) / 2, y = typeID+1, label = num), color = "black", size = 3)+
    #scale_y_continuous(breaks = 1:max(to_expand$typeID), labels = unique(CRE_granges$type)) +
    scale_fill_manual(values = cbp1, name = "CREs", labels = unique(CRE_granges$type)) +
    labs(x = paste0(chromosome, " position (bp)"), y = "CREs") +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),    # Remove y-axis ticks
      axis.text.y = element_blank(),     # Remove y-axis text
      axis.line.y = element_blank(),      # Remove y-axis line
      axis.ticks.x = element_blank(),    # Remove x-axis ticks
      axis.text.x = element_blank(),     # Remove x-axis text
      axis.line.x = element_blank(),      # Remove x-axis line
      axis.title.x = element_blank(),
    )
}


make_peaks_into_tracks <- function(annot_peaks, xrange){
  xmin <- xrange[1]
  xmax <- xrange[2]
  to_expand <- data.frame(ID = 1:nrow(annot_peaks),
                          start = annot_peaks$start,
                          end = annot_peaks$end)
  chromosome <- as.character(annot_peaks$seqnames[1])
  ggplot(to_expand, aes(xmin = start, xmax = end, ymin = 0, ymax =  0.8)) +
    geom_rect(color = "black", fill = "darkgray") + xlim(xmin, xmax) + ylim(0, 0.9)+
    geom_text(data = to_expand, aes(x = (start+end)/2, y = 0.4, label = ID), color = "black", size = 3)+
    labs(x = paste0(chromosome, " position (bp)"), y = "Peaks") +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),    # Remove y-axis ticks
      axis.text.y = element_blank(),     # Remove y-axis text
      axis.line.y = element_blank()      # Remove y-axis line
    )
}

make_genes_into_tracks <- function(annot_peaks, xrange){
  xmin <- xrange[1]
  xmax <- xrange[2]
  chromosome <- as.character(annot_peaks$seqnames[1])
  to_expand <- data.frame(start = max(annot_peaks$start, annot_peaks$geneStart),
                          end = min(annot_peaks$end, annot_peaks$geneEnd), gene_name = annot_peaks$SYMBOL)

  ggplot(to_expand, aes(xmin = start, xmax = end, ymin = 0, ymax =  0.8)) +
    geom_rect(color = "black", fill = "darkgray") + xlim(xmin, xmax) +
    geom_text(data = to_expand, aes(x = (start+end)/2, y = 0.4 , label = gene_name), color = "white", size = 3)+
    labs(x = paste0(chromosome, " position (bp)"), y = "Genes") +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),    # Remove y-axis ticks
      axis.text.y = element_blank(),     # Remove y-axis text
      axis.line.y = element_blank()      # Remove y-axis line
    )
}


make_transcript_into_tracks <- function(annot_peaks, xrange){
  xmin <- xrange[1]
  xmax <- xrange[2]
  chromosome <- as.character(annot_peaks$seqnames[1])
  transcript_info <- annot_peaks[, c("transcriptId", "SYMBOL", "geneStart", "geneEnd")]
  transcript_info <- transcript_info[!duplicated(transcript_info), ]
  to_expand <- data.frame(start = ifelse(xmin<transcript_info$geneStart, transcript_info$geneStart, xmin),
                          end = ifelse(xmax>transcript_info$geneEnd, transcript_info$geneEnd, xmax), transcript_name = transcript_info$transcriptId,
                          transcript_id =  1:nrow(transcript_info))

  ggplot(to_expand, aes(xmin = start, xmax = end, ymin =  transcript_id -0.4, ymax =   transcript_id +0.4)) +
    geom_rect(color = "black", fill = "white") + xlim(xmin, xmax) +
    geom_text(data = to_expand, aes(x = (start+end)/2, y =  transcript_id , label = transcript_name), color = "black", size = 2)+
    labs(x = paste0(chromosome, " position (bp)"), y = transcript_info$SYMBOL[1]) +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),    # Remove y-axis ticks
      axis.text.y = element_blank(),     # Remove y-axis text
      axis.line.y = element_blank()      # Remove y-axis line
    )
}

### work with more than two cell types
#' @importFrom patchwork wrap_plots
display_function <- function(all_peaks_summary, peaks.gr, annotation_using_ENCODE_CRE, chisq_test_on_CREs, iii){

  aa <- all_peaks_summary[[iii]]
  peak_info <- peaks.gr[iii]
  xmin <- start(peak_info)
  xmax <- end(peak_info)
  chromosome <- as.character(seqnames(peak_info))
  cut_summary <- melt(aa, id = "pos")
  colnames(cut_summary) <- c("pos", "celltype", "counts")

  p <- ggplot(
    data = cut_summary,
    mapping = aes(x = pos, y = counts, color = celltype)
  )
  p <- p + geom_point()  + facet_wrap(facets = ~celltype, strip.position = "left", ncol = 1, scales = "free")
  p <- p +
    scale_color_brewer(palette = "Set2") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("cleavage sites coverage"))
  p1 <- p + xlim(xmin, xmax) +  theme_classic()

  xrange <- c(xmin, xmax)

  p_CRE <- make_into_tracks(annotation_using_ENCODE_CRE[[iii]]$peaksENCODE, c(xmin, xmax))

  p_gene <- make_genes_into_tracks(annot_peaks[iii, ], c(xmin, xmax))

  freq <- chisq_test_on_CREs[[iii]]$merged
  freq <- apply(freq, 2, function(x){x/sum(x)})
  freq <- freq[-1, ]
  merged_fr <- data.frame(CRE_id = paste0("CRE", rownames(freq)), freq)
  merge_summary <- melt(merged_fr, id = "CRE_id")
  colnames(merge_summary)[2:3] <- c("celltype", "frequency")

  p_freq <- ggplot(merge_summary, aes(x = CRE_id, y = frequency, fill = celltype)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "CREs", y = "Frequency") +
    theme_classic() +
    scale_fill_brewer(palette = "Set2")  # Choosing a color palette

  patchwork::wrap_plots(p1, p_CRE, p_gene, p_freq, heights = c(6, 1, 0.5, 3), ncol = 1)

}




#' @importFrom viridis viridis_pal
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by filter
#' @importFrom patchwork wrap_plots plot_layout
display_a_peak <- function(for_plot_this_peak, celltypes = NULL){

  # Define a color palette for the cell types
  if(is.null(celltypes)){
    celltypes <- unique(for_plot_this_gene$cleavage$celltype)
  }

  if(length(celltypes) < 3) {
    colors <- viridis::viridis_pal(option = "D")(3)[1:length(celltypes)]
  } else {
    colors <- viridis::viridis_pal(option = "D")(length(celltypes))
  }

  celltype <- celltypes
  # Assign colors to cell types in a named vector
  color_mapping <- setNames(colors, celltype)

  anno <- for_plot_this_peak$annotation
  aa <- for_plot_this_peak$cleavage
  xmin <- anno$start
  xmax <- anno$end
  chromosome <- as.character(anno$seqnames)
  cut_summary <- melt(aa, id = "pos")
  colnames(cut_summary) <- c("pos", "celltype", "counts")

  if(!is.null(celltypes)){
    cut_summary <- cut_summary[cut_summary$celltype%in%celltypes, ]
  }

  cut_summary <- cut_summary[cut_summary$counts > 0, ]

  p <- ggplot(
    data = cut_summary,
    mapping = aes(x = pos, y = counts, color = celltype)
  )
  p <- p + geom_point()  + facet_wrap(facets = ~celltype, strip.position = "left", ncol = 1, scales = "free_y")
  p <- p +
    scale_color_manual(values = colors) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("cleavage sites coverage"))
  p1 <- p + xlim(xmin, xmax) +  theme_classic()

  xrange <- c(xmin, xmax)

  p_CRE <- make_into_tracks(for_plot_this_peak$CRE_info, c(xmin, xmax))

  p_gene <- make_genes_into_tracks(anno, c(xmin, xmax))

  freq <- for_plot_this_peak$freq_CRE

  p_freq <- ggplot(freq, aes(x = CRE_id, y = frequency, fill = celltype)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "CREs", y = "Frequency") +
    theme_classic() +
    scale_fill_manual(values = colors)

  patchwork::wrap_plots(p1, p_CRE, p_gene, p_freq, heights = c(6, 1, 0.5, 3), ncol = 1, guides = "collect") +
    patchwork::plot_layout(guides = "collect") +   # collect legends
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

}



#' @importFrom viridis viridis_pal
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by filter row_number n
#' @importFrom patchwork wrap_plots plot_layout
display_a_gene <- function(for_plot_this_gene, celltypes = NULL){

  # Define a color palette for the cell types
  if(is.null(celltypes)){
    celltypes <- unique(for_plot_this_gene$cleavage$celltype)
  }

  if(length(celltypes) < 3) {
    colors <- viridis::viridis_pal(option = "D")(3)[1:length(celltypes)]
  } else {
    colors <- viridis::viridis_pal(option = "D")(length(celltypes))
  }

  celltype <- celltypes
  color_mapping <- setNames(colors, celltype)


  anno <- for_plot_this_gene$peak_annotation
  aa <- for_plot_this_gene$cleavages
  xmin <- min(anno$start)
  xmax <- max(anno$end)
  chromosome <- anno$geneChr[1]
  cut_summary <- melt(aa, id = "pos")
  colnames(cut_summary) <- c("pos", "celltype", "counts")
  if(!is.null(celltypes)){
    cut_summary <- cut_summary[cut_summary$celltype%in%celltypes, ]
  }

  cut_summary <- cut_summary[cut_summary$counts > 0, ]

  p <- ggplot(
    data = cut_summary,
    mapping = aes(x = pos, y = counts, color = celltype)
  )
  p <- p + geom_point()  + facet_wrap(facets = ~celltype, strip.position = "left", ncol = 1, scales = "free_y")
  p <- p +
    scale_color_manual(values = color_mapping) +
    xlab(label = paste0("Chr", chromosome, " position (bp)")) +
    ylab(label = paste0("cleavage sites coverage"))
  p1 <- p + xlim(xmin, xmax) +  theme_classic()

  xrange <- c(xmin, xmax)

  p_CRE <- make_into_narrow_tracks(for_plot_this_gene$CRE_info, c(xmin, xmax))

  p_peak <- make_peaks_into_tracks(for_plot_this_gene$peak_annotation, c(xmin, xmax))

  p_transcript <- make_transcript_into_tracks(anno, c(xmin, xmax))

  CRE_freq <- for_plot_this_gene$CRE_freq_info
  CRE_freq$CRE_id <- as.numeric(gsub("CRE", "", CRE_freq$CRE_id))
  CRE_freq$peakID <- as.numeric(as.factor(CRE_freq$peakID))

  p_freq <- ggplot(CRE_freq, aes(x = CRE_id, y = frequency, fill = celltype)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "CREs", y = "Frequency")  +
    theme_classic() +
    scale_fill_manual(values = color_mapping)

  peak_CRE <- CRE_freq[, c(1, 2)]
  peak_CRE <- peak_CRE[!duplicated(peak_CRE),]
  colnames(peak_CRE) <- c("cre", "peak")
  extracted_df <- peak_CRE %>%
    dplyr::group_by(peak) %>%
    dplyr::filter(dplyr::row_number() %in% c(1, dplyr::n()))
  peak_num <- max(peak_CRE$peak)
  peak_CRE_draw <- data.frame(1:peak_num,
                              extracted_df[seq(1, peak_num*2, 2), 1],
                              extracted_df[seq(2, peak_num*2, 2), 1])
  colnames(peak_CRE_draw) <- c("peakID", "start", "end")

  p_freq <- p_freq + geom_segment(data = peak_CRE_draw,
                                  mapping = aes(x = start-0.25, xend = end + 0.25, y = 1, yend = 1), col = "black",
                                  inherit.aes = FALSE)
  p_freq <- p_freq + geom_segment(data = peak_CRE_draw, mapping = aes(x = start-0.52, xend = start-0.52, y = 0, yend = 1), col = "gray", inherit.aes = FALSE)

  p_freq <- p_freq + geom_text(data = peak_CRE_draw, mapping = aes(x = (start + end) / 2, y = 1, label = peakID), vjust = -0.2,  inherit.aes = FALSE)

  p_freq <- p_freq + scale_x_continuous(breaks = seq(1, max(CRE_freq$CRE_id))) +  # Customize x-axis ticks
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.5)) + geom_text(x =  0.1, y = 1.25, label = "Peak")

  num_CRE_type <- length(unique(for_plot_this_gene$CRE_info$type))
  num_transcripts <- length(unique(anno$transcriptId))
  num_celltypes <- length(celltypes)
  height_para <- c(2.5*num_celltypes, num_CRE_type, 1, 0.5*num_transcripts, 3)

  patchwork::wrap_plots(p1, p_CRE, p_peak,  p_transcript, p_freq, heights = height_para, ncol = 1, guides = "collect") +
  patchwork::plot_layout(guides = "collect") +   # collect legends
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

}



#' @importFrom viridis viridis_pal
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by filter
#' @importFrom patchwork wrap_plots
display_a_peak_with_coverage <- function(for_plot_this_peak, coverage_for_this_peak, celltypes = NULL){

  # Define a color palette for the cell types
  if(is.null(celltypes)){
    celltypes <- unique(for_plot_this_gene$cleavage$celltype)
  }

  if(length(celltypes) < 3) {
    colors <- viridis::viridis_pal(option = "D")(3)[1:length(celltypes)]
  } else {
    colors <- viridis::viridis_pal(option = "D")(length(celltypes))
  }

  celltype <- celltypes
  color_mapping <- setNames(colors, celltype)


  anno <- for_plot_this_peak$annotation
  aa <- for_plot_this_peak$cleavage
  xmin <- anno$start
  xmax <- anno$end
  chromosome <- as.character(anno$seqnames)
  cut_summary <- reshape2::melt(aa, id = "pos")
  colnames(cut_summary) <- c("pos", "celltype", "counts")
  if(!is.null(celltypes)){
    cut_summary <- cut_summary[cut_summary$celltype%in%celltypes, ]
  }

  cut_summary <- cut_summary[cut_summary$counts > 0, ]
  p <- ggplot(
    data = cut_summary,
    mapping = aes(x = pos, y = counts, color = celltype)
  )
  p <- p + geom_point()  + facet_wrap(facets = ~celltype, strip.position = "left", ncol = 1, scales = "free_y")
  p <- p +
    scale_color_manual(values = color_mapping) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("cleavage sites counts"))
  p1 <- p + xlim(xmin, xmax) +  theme_classic() + ylim(0, max(cut_summary$counts))

  xrange <- c(xmin, xmax)


  cut_summary <- reshape2::melt(coverage_for_this_peak, id = "pos")
  colnames(cut_summary) <- c("pos", "celltype", "counts")
  if(!is.null(celltypes)){
    cut_summary <- cut_summary[cut_summary$celltype%in%celltypes, ]
  }
  p <- ggplot(
    data = cut_summary,
    mapping = aes(x = pos, y = counts, color = celltype)
  )
  p <- p + geom_point()  + facet_wrap(facets = ~celltype, strip.position = "left", ncol = 1, scales = "free")
  p <- p +
    scale_color_manual(values = color_mapping) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("Coverage"))
  p2 <- p + xlim(xmin, xmax) +  theme_classic()
  p2

  p_CRE <- make_into_tracks(for_plot_this_peak$CRE_info, c(xmin, xmax))

  p_gene <- make_genes_into_tracks(anno, c(xmin, xmax))

  freq <- for_plot_this_peak$freq_CRE
  p_freq <- ggplot(freq, aes(x = CRE_id, y = frequency, fill = celltype)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "CREs", y = "Frequency") +
    theme_classic() +
    scale_fill_manual(values = color_mapping)  # Choosing a color palette

  num_CRE_type <- length(unique(for_plot_this_peak$CRE_info$type))
  num_celltypes <- length(celltypes)
  height_para <- c(2*num_celltypes, 2*num_celltypes, num_CRE_type,  0.5, 3)

  patchwork::wrap_plots(p2, p1, p_CRE, p_gene, p_freq, heights = height_para, ncol = 1, guides = "collect") +
      plot_layout(guides = "collect") +   # collect legends
      theme(plot.margin = margin(0, 0, 0, 0, "cm"))

}


