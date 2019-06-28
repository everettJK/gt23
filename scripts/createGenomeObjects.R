createRandomFragments <- function(refGenome, n = 1, fragWidth = 1){
  fragsPerChromosome <-  ceiling(n / length(refGenome))
  g <- unlist(GRangesList(lapply(1:length(refGenome), function(x){
    starts <- sample(1:(length(refGenome[[x]]) - fragWidth), fragsPerChromosome, replace = TRUE)
    GRanges(seqnames =  names(refGenome)[x], strand = '*', ranges = IRanges(start=starts, end = (starts + fragWidth - 1)))
  })))
  strand(g) <- sample(c('+', '-'), length(g), replace = TRUE)
  s <- as.character(BSgenome::getSeq(refGenome, g))
  s[! grepl('N', s)]
}

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm9)

hg38.randomFragments.100 <- createRandomFragments(refGenome = BSgenome.Hsapiens.UCSC.hg38, n = 500000, fragWidth = 100)
mm9.randomFragments.100  <- createRandomFragments(refGenome = BSgenome.Mmusculus.UCSC.mm9, n = 500000, fragWidth = 100)

save(hg38.randomFragments.100, file='../data/hg38.randomFragments.100.RData')
save(mm9.randomFragments.100, file='../data/mm9randomFragments.100.RData')