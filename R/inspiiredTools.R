addPositionID <- function(gr){
  gr$posid <- paste0(seqnames(gr), strand(gr), start(flank(gr, -1, start=T)))
  gr
}


#' Merge replicate samples and calculate abundances
#'
#' Merge each replicated sample {sampleName} into single samples {GTSP} and calculate clonal
#' abundances {estAbund}.
#'
#' @param gr A GRange object of intSites.
#' @param method Method ('fragments' or 'reads') to use to determine clonal abudnace.
#'
#' @return A GRange object of intSites where sequence ranges hanve been reduced to their start
#' positions, relicate samples have been combined into single samples and a clonal abundance
#' column has been added.
#'
#' @export
calcReplicateAbundances <- function(gr){

  if(!'posid' %in% names(mcols(gr))){
    gr <- addPositionID(gr)
  }

  if(!'reads' %in% names(mcols(gr))){
    stop("The provided GRange object does not have a 'reads' column")
  }

  if(!'GTSP' %in% names(mcols(gr))){
    stop("The provided GRange object does not have a 'GTSP' column")
  }

  names(gr) <- NULL
  gr$startPosition <- start(flank(gr, -1, start=T))

  d <- as.data.frame(gr)
  d$frags <- 0
  d$w <- paste0(d$sampleName, '/', d$width)

  o <- dplyr::group_by(d, GTSP, posid) %>%
    dplyr::mutate(reads = sum(reads),
                  frags = length(unique(w))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(estAbund = frags) %>%
    dplyr::group_by(GTSP) %>%
    dplyr::mutate(relAbund = (estAbund / sum(estAbund))*100) %>%
    dplyr::ungroup() %>%
    data.frame()

  o$start <- o$startPosition
  o$end   <- o$startPosition
  o$w <- NULL
  o$startPosition <- NULL
  
  GenomicRanges::makeGRangesFromDataFrame(o, keep.extra.columns=TRUE)
}


#' Standardize sequence ranges with Malani method
#'
#' @param gr A GRange object of intSites.
#' @param standardizeStart Boolean which determines if sequence range start positions should be standardized.
#' @param standardizeEnd Boolean which determines if sequence range end positions should be standardized.
#'
#' @return A GRange object of intSites where sequence ranges have been standardized.
#'
#' @export
standardizeSeqRanges.malani <- function(gr, standardizeStart = TRUE, standardizeEnd = TRUE, windowSize = 5L){

  # (!) clusterSites uses count() from matrixStats -- make sure that this version masks dplyr::count()

  if(standardizeStart){
    gr$Position <- ifelse(strand(gr) == "+", start(gr), end(gr))
    gr$Break    <- ifelse(strand(gr) == "+", end(gr),   start(gr))
    gr$Score    <- 95
    gr$qEnd     <- width(gr)

    gr.std <- hiReadsProcessor::clusterSites(psl.rd = gr, weight = rep(1, length(gr)), windowSize = windowSize)

    start(gr.std) <- ifelse(strand(gr.std) == "+", gr.std$clusteredPosition, gr.std$Break)
    end(gr.std)   <- ifelse(strand(gr.std) == "-", gr.std$clusteredPosition, gr.std$Break)

    for (i in c('Position','Break','Score', 'qEnd', 'clusteredPosition',
                'clonecount', 'clusterTopHit')){ mcols(gr.std[,i]) <- NULL }

    gr <- gr.std
  }

  if(standardizeEnd){
    gr$Position <- ifelse(strand(gr) == "-", start(gr), end(gr))
    gr$Break    <- ifelse(strand(gr) == "-", end(gr),   start(gr))
    gr$Score    <- 95
    gr$qEnd     <- width(gr)

    gr.std <- hiReadsProcessor::clusterSites(psl.rd = gr, weight = rep(1, length(gr)), windowSize = windowSize)

    start(gr.std) <- ifelse(strand(gr.std) == "-", gr.std$clusteredPosition, gr.std$Break)
    end(gr.std)   <- ifelse(strand(gr.std) == "+", gr.std$clusteredPosition, gr.std$Break)

    for (i in c('Position','Break','Score', 'qEnd', 'clusteredPosition',
                'clonecount', 'clusterTopHit')){ mcols(gr.std[,i]) <- NULL }
    gr <- gr.std
  }

  gr
}


