chunkingVector <- function(x, n) ceiling(seq_along(x)/(length(x)/n))



parse_cdhitest_output <- function(file){
  # cd-hit-est must be ran with -d 0 flag in order to preserve seq names up to first white space.
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, '>Cluster'))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x){
     gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(x, '>([^\\s]+)'), function(y){ y[,2] })))
  })
}


parsePairedEndIndexedSeqRun <- function(sampleFile, I1file, R1file, R2file, CPUs=1, 
                                        I1trimFun=NULL, R1trimFun=NULL, R2trimFun=NULL){
  
  sampleData <- read.csv(sampleFile, stringsAsFactors = FALSE)
  
  message('Reading I1 ...'); I1 <- ShortRead::readFastq(I1file)
  message('Reading R1 ...'); R1 <- ShortRead::readFastq(R1file)
  message('Reading R2 ...'); R2 <- ShortRead::readFastq(R2file)
  
  if(! is.null(I1trimFun)) I1 <- I1trimFun(I1)
  if(! is.null(R1trimFun)) R1 <- I1trimFun(R1)
  if(! is.null(R2trimFun)) R2 <- I1trimFun(R2)
  
  R1.df <- data.frame(sequence = as.character(ShortRead::sread(R1)))
  R2.df <- data.frame(sequence = as.character(ShortRead::sread(R2)))
  
  samples <- sampleData[match(as.character(ShortRead::sread(I1)), sampleData$bcSeq),]$alias
  samples[is.na(samples)] <- 'unknown'
  
  R1.lst <- split(R1.df, samples)
  R2.lst <- split(R2.df, samples)
  
  R1.df <- plyr::ldply(R1.lst, .id='sampleName')
  R2.df <- plyr::ldply(R2.lst, .id='sampleName')
  
  R1.df$read <- 'R1'
  R2.df$read <- 'R2'
  
  d <- dplyr::bind_rows(R1.df, R2.df)
  
  message('Tallying unique reads ...')
  cluster <- parallel::makeCluster(CPUs)
  o <- dplyr::bind_rows(parallel::parLapply(cluster, split(d, d$sampleName), function(y){
    library(dplyr)
    dplyr::group_by(y, sequence) %>%
      dplyr::summarise(sampleName = sampleName[1], nReads = n()) %>%
      dplyr::ungroup() %>%
      data.frame()
  }))
  parallel::stopCluster(cluster)
  
  o
}


getSiteFlankingSequences <- function(g, n=12){
  library(GenomicRanges)
  g.upstream   <- flank(g, n, start = TRUE)
  g.downstream <- flank(g, n, start = FALSE)
  
  # flank() does not include the genomic range that its flanking.
  # Here we shift the down stream ranges 1 NT so that the upstream and downstream ranges are contiguous.
  i <- which(strand(g.downstream) == '-')
  g.downstream[i]  <- shift(g.downstream[i], 1)
  g.downstream[-i] <- shift(g.downstream[-i], -1)
  
  g.upstream.seq   <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, g.upstream)
  g.downstream.seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, g.downstream)
  
  # Name the sequences for future MSA analysis.
  names(g.upstream.seq) <- paste0('s', 1:length(g.upstream.seq))
  names(g.downstream.seq) <- paste0('s', 1:length(g.downstream.seq))
  
  list('upstream' = g.upstream.seq, 'downstream' = g.downstream.seq)
}


# Function from online source.
# https://stackoverflow.com/questions/11340444
numShortHand <- function (number, rounding=F, digits=ifelse(rounding, NA, 2))
{
  lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
           0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
           1e+24, 1e+27)
  pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k",
           "M", "G", "T", "P", "E", "Z", "Y", NA)
  ix <- findInterval(number, lut)
  if (ix>0 && ix<length(lut) && lut[ix]!=1) {
    if (rounding==T && !is.numeric(digits)) {
      sistring <- paste(round(number/lut[ix]), pre[ix])
    }
    else if (rounding == T || is.numeric(digits)) {
      sistring <- paste(signif(number/lut[ix], digits), pre[ix])
    }
    else {
      sistring <- paste(number/lut[ix], pre[ix])
    }
  }
  else {
    sistring <- as.character(number)
  }
  
  sistring <- gsub('\\s', '', sistring)
  return(sistring)
}


englishList <- function(v, conjunction = 'and'){
  if(length(v) == 1) return(v)
  if(length(v) == 2) return(paste(v, collapse = paste0(' ', conjunction, ' '))) 
  paste0(paste0(v[1:length(v)-1], collapse = ', '), paste0(' ', conjunction, ' '), v[length(v)])
}



ppNum <- function(n) format(n,big.mark=",",scientific=FALSE, trim=TRUE)



createColorPalette <- function(n){
  library(RColorBrewer)
  library(grDevices)
  colorRampPalette(brewer.pal(12, "Paired"))(n)
}


# Create 30 alphanumeric character tmp file name.
tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }


# Determine the frequency of strings in a vector of strings and return a 
# numeric vector ordering the strings by decreasing frequency.
#   stringFreqOrder(c('b', 'a', 'a', 'c', 'a', 'b'))
#    2 1 1 3 1 2
stringFreqOrder <- function(s){
  i <- names(sort(table(s), decreasing = TRUE))
  match(s, i)
}



# Sum a numeric vector ignoring NAs - return NA if all values are NA.
NA.sum <- function(x){
  if(all(is.na(x))){
    return(NA)
  } else {
    return(sum(x, na.rm = TRUE))
  }
}




# Helper function to wait for a file to appear -- helpful for waiting for system calls to create outputs.
waitForFile <- function(f, seconds = 1){
  repeat
  {
    if(file.exists(f)) break
    Sys.sleep(seconds)
  }
  
  return(TRUE)
}



# Create a string denoting the occurances of strings in a vector of strings.
# {n} denoted the tally value after which no more is reported.
# tallyString(c('a', 'a', 'a', 'b', 'a', 'b', 'c'), n = 2)
# 'a x4, b x2 +1 more'

tallyString <- function(v, n){ 
  v <- sort(table(v), decreasing = TRUE); 
  if(length(v) <= n){
    n <- length(v) 
    additional <- ''
  } else {
    additional <- paste0(' +', length(v) - n, ' more')
  }
  
  v2 <- v[1:n]
  paste0(paste0(paste0(names(v2), ' x', v2), collapse = ', '), additional) 
}




