#' Standardize genomic fragments from intSite experiments using the gintools package.
#'
#' @param frags GenomicRange object containing genomic fragments where fragments with positive strands 
#' are assumed to have their integration positions on the left and break points on the right (opposite assumed from negative strand).
#'
#' @param CPUs Number of CPUs to use for break point standardizations.
#' @param countsCol Name of meta data column containing the number of times fragments were observed. 
#'
#' @importFrom magrittr '%>%'
#' @export
stdIntSiteFragments <- function(frags, CPUs = 10, countsCol = 'reads'){ 
  # Setup parallelization.
  cluster <- parallel::makeCluster(CPUs)
  parallel::clusterExport(cl = cluster, envir = environment(), varlist = c('countsCol'))

  # Standardize fragment start positions.
  frags <- gintools::standardize_sites(frags, counts.col = countsCol)


   # Standardize break points by sample.
   frags <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, as.list(split(frags, frags$sampleName)), function(x){
              gintools::refine_breakpoints(x, counts.col = countsCol)
            })))

   # Stop the cluster.
   parallel::stopCluster(cluster)

   # Create a postion id now that fragment boundaries have been standardized.
   frags$posid <- paste0(GenomicRanges::seqnames(frags), GenomicRanges::strand(frags), GenomicRanges::start(GenomicRanges::flank(frags, -1, start = T)))


   # Now that the fragment positions have been adjusted, merge ranges and re-tally the counts.
   dplyr::group_by(data.frame(frags), cellType, timePoint, start, end, strand) %>%
   dplyr::mutate(reads = sum(reads)) %>%
   dplyr::slice(1) %>%
   dplyr::ungroup() %>%
   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}



#' Collapse technical replicates and estimate abundances by number of unique fragments associated with integration positions.
#'
#' @param f GRange object of standardized genomic fragments.
#' @return GRange object of integration positions including estimated abundances.
#' @importFrom magrittr '%>%'
#' @export
collapseReplicatesCalcAbunds <- function(f){
  # Conversion of GRange object to data frame creates width column.
  f <- data.frame(f)
  f$start <- ifelse(as.character(f$strand) == '+', f$start, f$end)
  f$end   <- f$start
  dplyr::group_by(f, GTSP, posid) %>%
  dplyr::mutate(reads = sum(reads)) %>%
  dplyr::mutate(estAbund = dplyr::n_distinct(width)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sampleName) %>%
  dplyr::group_by(GTSP) %>%
  dplyr::mutate(relAbund = (estAbund / sum(estAbund))*100) %>%
  dplyr::ungroup() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}
  
 


#' Create a list of default file mappings.
#' 
#' An empty character vector object named 'none' is available for instances where 
#' gene lists are not available.
#'
#' @return List of default gt23 file mappings.
#'
#' @export
defaultGenomeFileMappings <- function(){
   list('hg38'    = list('genes'             = 'hg38.refSeqGenesGRanges', 
                         'exons'             = 'hg38.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'),
        'hg18'    = list('genes'             = 'hg18.refSeqGenesGRanges', 
                         'exons'             = 'hg18.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                        'lymphomaGenesList'  = 'hg38.lymphomaGenesList'),
        'mm9'     = list('genes'             = 'mm9.refSeqGenesGRanges', 
                         'exons'             = 'mm9.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'mm9.oncoGeneList',
                         'lymphomaGenesList' = 'none'),
        'canFam3' = list('genes'             = 'canFam3.humanXeno.refSeqGenesGRanges', 
                         'exons'             = 'canFam3.humanXeno.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'),
        'macFas5' = list('genes'             = 'macFas5.humanXeno.refSeqGenesGRanges', 
                         'exons'             = 'macFas5.humanXeno.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'))
        
 }
 
 
#' Validate that each named data object in genome file mapping list is present in the gt23 package.
#'
#' The special object 'none' is allowed which is an empty character vector object.
#'
#' @return Error if one or more named data objects are not present.
#'
#' @export
 validateGenomeFileMap <- function(x){
   availableObjects <- data(package='gt23')$results[,3]
   
   invisible(lapply(x, function(i){
     i <- unlist(i)
     i <- i[! is.na(i)]
     if(any(! i %in% availableObjects)) 
       stop(paste0(paste0(i[! i %in% availableObjects], collapse = ', '), ' could not be found in the package database. ',
                   'Availabe data objects include: ', paste0(availableObjects, collapse = ', ')))
   }))
 }
 
 
 #' Calculate nearest genomic and onocogene features. 
 #'
 #' @return GRange object updated with genomic feature annotations.
 #'
 #' @export
 annotateIntSites <- function(sites, CPUs = 20, genomeFileMap = NULL){
   
   if(is.null(genomeFileMap)) genomeFileMap <- gt23::defaultGenomeFileMappings()
   gt23::validateGenomeFileMap(genomeFileMap)
   
   # Here we convert GRanges to data frames after annoations are added and then back to a single 
   # GRange object so that different metadata cols can be concatenated where nonexistinting columns will be NA.
   
   GenomicRanges::makeGRangesFromDataFrame(dplyr::bind_rows(lapply(split(sites, paste(sites$patient, sites$refGenome)), function(x){
     
     message(paste('Starting', x$patient[1], '/', x$refGenome[1]))
     
     # Use the refGenome id to retrieve the appropriate gene boundary data.
     if(! x$refGenome[1] %in% names(genomeFileMap)) stop(paste('Error -- ', x$refGenome[1], ' is not defined in the genomeFileMap list.'))
                                                         
     genome_refSeq      <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$genes)))
     genome_refSeqExons <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$exons)))
     oncoGeneList       <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$oncoGeneList)))
     lymphomaGenesList  <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$lymphomaGenesList)))
 
     # Setup parallelization.
     cluster <- parallel::makeCluster(CPUs)
     parallel::clusterExport(cl=cluster, envir = environment(), varlist = c('CPUs', 'genome_refSeq', 'genome_refSeqExons', 'oncoGeneList', 'lymphomaGenesList'))
    
     # Create a splitting vector for parallelization.
     x$s <- dplyr::ntile(seq_along(x), CPUs)
     names(x) <- NULL
     
     x <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(x, x$s), function(x2){
          library(GenomicRanges)
          library(gt23)
          library(dplyr)
       
          # Create an order column to ensure that the incoming ranges keep the same order.
          x2$n <- 1:length(x2)
       
          # Nearest gene boundary
          x2 <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons)
          x2 <- x2[order(x2$n)]
       
         # Nearest oncogene
         if(length(oncoGeneList) > 0){
           o <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons, geneList=gt23::humanOncoGenesList)
           o <- o[order(o$n)]
       
           d <- data.frame(GenomicRanges::mcols(x2))
           d <- d[order(d$n),]
           stopifnot(all(o$n == d$n))
           d$nearestOncoFeature       <- o$nearestFeature
           d$nearestOncoFeatureDist   <- o$nearestFeatureDist
           d$nearestOncoFeatureStrand <- o$nearestFeatureStrand
           GenomicRanges::mcols(x2) <- d
         }
     
         # Nearest lymphoma gene
         if(length(lymphomaGenesList) > 0){
           o <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons, geneList=gt23::humanLymphomaGenesList)
           o <- o[order(o$n)]
       
           d <- data.frame(GenomicRanges::mcols(x2))
           d <- d[order(d$n),]
           stopifnot(all(o$n == d$n))
           d$nearestlymphomaFeature       <- o$nearestFeature
           d$nearestlymphomaFeatureDist   <- o$nearestFeatureDist
           d$nearestlymphomaFeatureStrand <- o$nearestFeatureStrand
           GenomicRanges::mcols(x2) <- d
         }
       
         x2$n <- NULL
         x2
     })))
     
     x$s      <- NULL
     names(x) <- NULL
     parallel::stopCluster(cluster)
     data.frame(x)
   })), keep.extra.columns = TRUE)
 }
 
 
 
 
 

