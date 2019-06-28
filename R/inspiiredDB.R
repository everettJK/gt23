#' Retrieve patient sequencing run data
#'
#' Create a dataframe of sequencing run ids associated with 1 or more patient identifiers.
#'
#' @param patients A vector of patient identifiers.
#' @param specimenManagementDBgroup MySQL database group for specimen management database.
#' @param intsitesDBgroup MySQL database group for intSite database.
#'
#' @return A data frame with both run ids and run dates.
#'
#' @export
getPatientRunDetails <- function(patients, specimenManagementDBgroup='specimen_management', intsitesDBgroup='intsites_miseq'){
  
  # Create DB connections
  dbConn1 <- DBI::dbConnect(MySQL(fetch.default.rec = Inf), group=specimenManagementDBgroup)
  dbConn2 <- DBI::dbConnect(MySQL(fetch.default.rec = Inf), group=intsitesDBgroup)

  # Retrieve GTSP ids
  query <- paste('select SpecimenAccNum from gtsp where Patient in (', paste(sQuote(patients), collapse=', '), ')')
  GTSPs <- unique(unname(unlist(DBI::dbGetQuery(dbConn1, query))))

  # Retrieve run ids
  query <- paste('select miseqid from samples where sampleName like ', paste(sQuote(paste0(GTSPs, '-%')), collapse=' or sampleName like '))
  runIDs <- unique(unname(unlist(DBI::dbGetQuery(dbConn2, query))))

  DBI::dbDisconnect(dbConn1)
  DBI::dbDisconnect(dbConn2)

  d <- data.frame(runID=runIDs,
                  runDate=as.Date(str_extract_all(runIDs, '(^[^_]+)', simplify = TRUE), "%y%m%d"))

  d[order(d$runDate),]
}





#' Retrieve sequencing run patient details
#'
#' Create a dataframe of detailing a specific sequencing run.
#'
#' @param runID A single MiSeq sequencing run id.
#' @param specimenManagementDBgroup MySQL database group for specimen management database.
#' @param intsitesDBgroup MySQL database group for intSite database.
#'
#' @return A data frame detailing patients and sample details.
#'
#' @export
getRunDetails <- function(runID, specimenManagementDBgroup='specimen_management', intsitesDBgroup='intsites_miseq'){

  # Create DB connections
  dbConn1 <- DBI::dbConnect(MySQL(fetch.default.rec = Inf), group=specimenManagementDBgroup)
  dbConn2 <- DBI::dbConnect(MySQL(fetch.default.rec = Inf), group=intsitesDBgroup)

  sampleDB <- DBI::dbGetQuery(dbConn1, 'select * from gtsp')

  # Retrieve run ids
  query <- paste0('select distinct sampleName from samples where miseqid="', runID, '"')
  sampleNames <- unname(unlist(DBI::dbGetQuery(dbConn2, query)))
  sampleNames <- unique(sub('\\-\\d+$', '', sampleNames[grep('^GTSP', sampleNames)]))

  DBI::dbDisconnect(dbConn1)
  DBI::dbDisconnect(dbConn2)

  sampleData <- data.frame(GTSP=sampleNames)
  sampleData$patient <- sampleDB[match(sampleData$GTSP, sampleDB$SpecimenAccNum),]$Patient
  sampleData$trial <- sampleDB[match(sampleData$GTSP, sampleDB$SpecimenAccNum),]$Trial
  sampleData$timePoint <- sampleDB[match(sampleData$GTSP, sampleDB$SpecimenAccNum),]$Timepoint
  sampleData$cellType <- sampleDB[match(sampleData$GTSP, sampleDB$SpecimenAccNum),]$CellType
  sampleData
}



#' Disconnect all MySQL database connections.
#'
#' @return nothings
#' 
#' @export
disconnectAllDBs <- function(){
  invisible(lapply(DBI::dbListConnections(RMySQL::MySQL()), DBI::dbDisconnect))
}


#' Parse the sample replicate ids in a sequencing sample sheet.
#'
#' @return Data frame of unique samples including patient, cell type and time point information.
#' 
#' @export
getSampleSheetDetails <- function(sampleSheetFile, specimenManagementDBgroup='specimen_management'){
  dbConn1  <- DBI::dbConnect(RMySQL::MySQL(), group=specimenManagementDBgroup)
  sampleDB <- DBI::dbGetQuery(dbConn1, 'select * from gtsp')
  DBI::dbDisconnect(dbConn1)
  f <- read.table(file = sampleSheetFile, sep = ',', header = TRUE)
  s <- sub('\\-\\d+$', '', as.character(f$alias))
  s <- unique(s[grepl('GTSP', s)])
  dplyr::filter(sampleDB, SpecimenAccNum %in% s) %>%
  dplyr::select(Trial, Patient, CellType, Timepoint) %>%
  dplyr::distinct()
}

  

#' Retreive intSite data
#'
#' Retreive inSite data for either a number of patients, sequencing runs or samples.
#'
#' @param sampleDB.group Sample database connection group
#' @param intSiteDB.group IntSite database connection group
#' @param patients Vector of patient IDs for which to retrieve intSites.
#' @param samples Vector of samples IDs for which to retrieve intSites.
#' @param runIDs Vector of sequencing run ids for which to retrieve intSites.
#' @param roundTimePoints Logical determining if time points should be rounded.
#'
#' @return A GRange object of intSites.
#'
#' @export
getIntSiteData <- function(samples, sampleDB.group, intSiteDB.group,
                           roundTimePoints=TRUE, runIDs=NULL){
  options(useFancyQuotes = FALSE)

  if(length(samples) == 0) stop('One or more sample ids must be provided.')
  
  dbConn1  <- DBI::dbConnect(RMySQL::MySQL(), group=sampleDB.group)
  dbConn2  <- DBI::dbConnect(RMySQL::MySQL(), group=intSiteDB.group)
  
  intSiteSamples <- DBI::dbGetQuery(dbConn2, 'select * from samples')
  intSiteSamples$GTSP <- gsub('\\-\\d+$', '', intSiteSamples$sampleName)
  
  # Retreive intSite db numeric sample ids.
  sampleIDs <- unique(base::subset(intSiteSamples, GTSP %in% samples)$sampleID)
  
  if(length(sampleIDs) == 0) stop('Error: no intSite sample ids have been selected.')
  
  replicateQuery <- paste('samples.sampleID in (', paste0(sampleIDs, collapse = ','), ')')
    
  if(!is.null(runIDs)){
    runIDs <- paste0('(miseqid in (', paste0(sQuote(runIDs), collapse = ','), ')) and ')
  } else {
    runIDs = ''
  }
  
  q <- sprintf("select position, chr, strand, breakpoint, count, refGenome,
               sampleName from sites left join samples on
               sites.sampleID = samples.sampleID
               left join pcrbreakpoints on
               pcrbreakpoints.siteID = sites.siteID
               where %s (%s)", runIDs, replicateQuery)
  
  options(warn=-1)
  sampleData <- DBI::dbGetQuery(dbConn1, "select * from gtsp")
  sites <- dbGetQuery(dbConn2, q)
  options(warn=0)
  
  if(nrow(sites) == 0) return(GRanges())
  
  sites$GTSP      <- as.character(sub('\\-\\d+', '', sites$sampleName))
  sites$patient   <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'Patient']
  sites$timePoint <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'Timepoint']
  sites$cellType  <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'CellType']
  
  sites$timePoint <- toupper(sites$timePoint)
  sites$timePoint <- gsub('_', '.', sites$timePoint)
  
  sites$timePointType<- stringr::str_match(sites$timePoint, '[DMY]')
  sites$timePointType[which(is.na(sites$timePointType))] <- 'X'
  
  sites <- dplyr::bind_rows(lapply(split(sites, sites$timePointType), function(x){
    n <- as.numeric(stringr::str_match(x$timePoint, '[\\d\\.]+')) * ifelse(grepl('\\-', x$timePoint), -1, 1)
    
    if(x$timePointType[1] == 'D'){
      x$timePointMonths <- n / 30.4167
      x$timePointDays   <- n
    } else if(x$timePointType[1] == 'M'){
      x$timePointMonths <- n
      x$timePointDays   <- n * 30.4167
    } else if(x$timePointType[1] == 'Y'){
      x$timePointMonths <- n * 12
      x$timePointDays   <- n * 365
    } else {
      message('Warning - could not determine date unit for: ', x$timePointType[1])
      x$timePointMonths <- n
      x$timePointDays   <- n 
    }
    x
  }))
  sites$timePointType <- NULL
  
  if(roundTimePoints) sites$timePointMonths <- base::round(sites$timePointMonths)
  
  intSites <-GenomicRanges::GRanges(seqnames = S4Vectors::Rle(sites$chr),
                         ranges = IRanges::IRanges(start = pmin(sites$position, sites$breakpoint),
                                                   end   = pmax(sites$position, sites$breakpoint)),
                         strand = S4Vectors::Rle(sites$strand))
  
  sites <- dplyr::rename(sites, reads = count)
  mcols(intSites) <- data.frame(dplyr::select(sites, refGenome, reads, patient, sampleName, GTSP, 
                                              cellType, timePoint, timePointDays, timePointMonths))
  intSites
}
