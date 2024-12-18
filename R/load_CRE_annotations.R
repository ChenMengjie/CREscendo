load_CRE_annotations <- function(genome = "hg38") {
  if(genome%in%"hg38"){
    file_path <- system.file("extdata", "hg38_cCREs.rds", package = "CREscendo")
  }
  if(genome%in%"mm10"){
    file_path <- system.file("extdata", "mm10_cCREs.rds", package = "CREscendo")
  }
  if (file_path == "") stop("Annotation file not found in the package.")
  gr <- readRDS(file_path)
  return(gr)
}
