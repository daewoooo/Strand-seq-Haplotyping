#' This function localize seqments of continous values (zeros and ones)
#' 
#' @importFrom fastseg fastseg
#'
#' @author David Porubsky
 
findRecombs <- function(comparisons=NULL, minSeg=0, smooth=3, chromosome=NULL, filtSize=5000000) {

  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1  
    } 
  }
  
  collapseBins <- function(gr, id.field=0) {
    ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
  }
  
  #localize segments in 0/1 vector
  segments <- GRangesList()
  for (i in 1:length(comparisons)) {
    comparison <- comparisons[[i]]
    index <- names(comparisons[i])
    comp.vector <- as.numeric(comparison[,2])
    
    if (length(comp.vector) > 2*minSeg) {
      segs <- fastseg(comp.vector, minSeg=minSeg)
    
      while (any(segs$num.mark <= smooth)) {
        toSwitch <- which(segs$num.mark <= smooth)
        switch.segs <- segs[toSwitch]
        switch.pos <- mapply(function(x,y) {x:y}, x=switch.segs$startRow, y=switch.segs$endRow)
        switch.pos <- unlist(switch.pos)
    
        switched.vals <- sapply(comp.vector[switch.pos], switchValue) #SWITCH
        comp.vector[switch.pos] <- switched.vals
        segs <- fastseg(comp.vector, minSeg=minSeg)
      }
      gen.ranges <- IRanges(start=comparison[,1][segs$startRow], end=comparison[,1][segs$endRow])
      ranges(segs) <- gen.ranges
      segs$index <- index
      
      segs$match[segs$seg.mean <= 0.25] <- 'hap1'
      segs$match[segs$seg.mean >= 0.75] <- 'hap2'
      segs$match[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'mix'
      
      segments[[index]] <- segs
    }
  }	
  
  recombs <- GRangesList()
  for (i in 1:length(segments)) {
    segm <- segments[[i]]
    index <- names(segments[i])
    
    if (length(segm) > 1) {
      segm <- collapseBins(gr = segm, id.field = 7)
      segm <- segm[width(segm) >= filtSize]
      segm <- collapseBins(gr = segm, id.field = 7)	
      segments[[i]] <- segm
    }	 
    
    if (length(segm) > 1) {
      suppressWarnings( recomb <- gaps(segments[[i]], start = start(segments[[1]])) )
      index <- names(segments[i])
      start(recomb) <- start(recomb)-1
      end(recomb) <- end(recomb)+1
      recomb$index <- index
      recombs[[index]] <- recomb
    }  
  }
  
  seqments.gr <- unlist(segments)
  names(seqments.gr) <- NULL
  segments.df <- as(seqments.gr, "data.frame")
  segments.df <- data.frame(chromosome=chromosome, start=segments.df$start, end=segments.df$end, seg.mean=segments.df$seg.mean, index=segments.df$index, match=segments.df$match)
  
  recombs.gr <- unlist(recombs)
  names(recombs.gr) <- NULL
  recombs.df <- as(recombs.gr, "data.frame")
  recombs.df <- data.frame(chromosome=chromosome, start=recombs.df$start, end=recombs.df$end, range=recombs.df$width, index=recombs.df$index)
  
  results <- list()
  results[['segments']] <- segments.df
  results[['recombs']] <- recombs.df
  return(results)
}


################################################################################################################################################################################


#' This function localize seqments of continous values (zeros and ones)
#' 
#' @param infile File containing haplotype comparisons between one parent and a child (see example below).
#' @param outputDirectory Folder to output the results. If it does not exist it will be created.
#' @param minSeg Length of a the minimal segment to consider
#' @param smooth 
#' @param filtSize Filter out segments smaller than 
#' @param chromosome
#'
#' @author David Porubsky

#Example infile
#P1 - paternal allele 1
#P2 - paternal allele 2
#C1 - child allele 1
#C2 - child allele 2
#C1vsP1/P2 - comparison of C1 with P1 and P2
#C2vsP1/P2 - comparison of C2 with P1 and P2

#Chr	Pos	P1	P2	C1	C2	C1vsP1/P2	C2vsP1/P2
#1	695745	G	A	G	A	P1	P2
#1	787188	C	T	T	T	P2	P2
#1	791524	T	C	C	-	P2	NA
#1	791543	C	A	A	A	P2	P2
#1	791721	T	C	C	C	P2	P2
#1	791806	C	T	T	-	P2	NA

localizeRECs <- function(infile, outputDirectory="./Recomb_analysis", minSeg=100, smooth=2, filtSize = 5000000, chromosome=NULL) {
  
  if (!file.exists(outputDirectory)) {
    dir.create(outputDirectory)  
  }
  
  results <- file.path(outputDirectory, infile)
  if (!file.exists(results)) {
    dir.create(results)  
  }
  
  #read in data and split homologue comparisons to separate vectors
  data <- read.table(infile, header=T, stringsAsFactors = F)
  C1vsP1.P2 <- data[,c('Pos', 'C1vsP1.P2')]
  C2vsP1.P2 <- data[,c('Pos', 'C2vsP1.P2')]
  C1vsP1.P2 <- C1vsP1.P2[!is.na(C1vsP1.P2['C1vsP1.P2']) & C1vsP1.P2['C1vsP1.P2'] != 'E',]
  C2vsP1.P2 <- C2vsP1.P2[!is.na(C2vsP1.P2['C2vsP1.P2']) & C2vsP1.P2['C2vsP1.P2'] != 'E',]
  
  #encode every comparison to 0's and 1's
  C1vsP1.P2$C1vsP1.P2[C1vsP1.P2$C1vsP1.P2 == 'F1'] <- 0
  C1vsP1.P2$C1vsP1.P2[C1vsP1.P2$C1vsP1.P2 == 'F2'] <- 1
  C2vsP1.P2$C2vsP1.P2[C2vsP1.P2$C2vsP1.P2 == 'F1'] <- 0
  C2vsP1.P2$C2vsP1.P2[C2vsP1.P2$C2vsP1.P2 == 'F2'] <- 1
  
  comparisons <- list()
  comparisons[['C1vsP1.P2']] <- C1vsP1.P2
  comparisons[['C2vsP1.P2']] <- C2vsP1.P2
  
  rec.segm <- findRecombs(comparisons=comparisons, minSeg=minSeg, smooth=smooth, chromosome=chromosome, filtSize=filtSize)
  
  segments <- rec.segm$segments
  segments <- split(segments, segments$index)
  
  for (i in 1:length(segments)) {
    write.table(segments[[i]], file=file.path(results, paste0(names(segments[i]), "_segments.txt")), row.names=F)
  }
  
  recombs <- rec.segm$recombs
  recombs <- split(recombs, recombs$index)
  
  for (i in 1:length(recombs)) {
    write.table(recombs[[i]], file=file.path(results, paste0(names(recombs[i]), "_breaks.txt")), row.names=F)
  }
}


#TO RUN THE FUNCTIONS
#load the functions (findRecombs and localizeRECs) into the R environment
#then run this command
localizeRECs(<infile>, outputDirectory="./Recomb_analysis", minSeg=100, smooth=2, filtSize = 5000000, chromosome="1")

