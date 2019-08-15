library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm9)
set.seed(46)

createRandomFragments <- function(refGenome, n = 1, fragWidth = 1){
  
  min_seqs <- grep("_", seqnames(refGenome), fixed=TRUE, invert=TRUE, value=TRUE)
  refGenome@user_seqnames <- setNames(min_seqs, min_seqs)
  refGenome@seqinfo <- refGenome@seqinfo[min_seqs]
  
  fragsPerChromosome <-  ceiling(n / length(refGenome))
  g <- unlist(GRangesList(lapply(1:length(refGenome), function(x){
    starts <- sample(1:(length(refGenome[[x]]) - fragWidth), fragsPerChromosome, replace = TRUE)
    GRanges(seqnames =  names(refGenome)[x], strand = '*', ranges = IRanges(start=starts, end = (starts + fragWidth - 1)))
  })))
  
  p <- paste0(seqnames(g), ':', start(g), '-', end(g))
  s <- as.character(BSgenome::getSeq(refGenome, g))
  i <- grepl('N', s)
  p <- p[!i]
  s <- s[!i]
  names(s) <- p
  s
}

hg38.randomFragments.1000000.100 <- sample(createRandomFragments(refGenome = BSgenome.Hsapiens.UCSC.hg38, n = 1500000, fragWidth = 100), 1000000)
mm9.randomFragments.1000000.100  <- sample(createRandomFragments(refGenome = BSgenome.Mmusculus.UCSC.mm9, n = 1500000, fragWidth = 100), 1000000)

save(hg38.randomFragments.1000000.100, file='../data/hg38.randomFragments.1000000.100.RData', compress = TRUE, compression_level = 9)
save(mm9.randomFragments.1000000.100, file='../data/mm9randomFragments.1000000.100.RData', compress = TRUE, compression_level = 9)