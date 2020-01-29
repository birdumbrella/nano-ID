

### prewd ###


prewd = file.path(getwd(),"...")


### libraries ###


library(Biostrings)
library(foreach)
library(doParallel)


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects","fast5.file","fast5.files")


### load objects ###


load(file.path("...","fast5.files.RData"))


### alignment.reconstruction.function ###


alignment.reconstruction.function = function(read.name){
  fast5.name = read.to.fast5.name[read.name]
  read = read.sequence.list[[read.name]]
  which.strand = as.character(strand(bam[read.name]))
  if (which.strand == "-"){read = as.character(reverseComplement(as(read,'DNAString')))}
  reference = aligned.sequence.list[[read.name]]
  cigar = cigar.list[[read.name]]
  exploded.cigar = as.vector(Rle(names(cigar),cigar))
  alignment.length = length(exploded.cigar)
  aligned.read = rep("-",alignment.length)
  aligned.read[which(exploded.cigar %in% c("M","I","S","H"))] = strsplit(read,split = "")[[1]]
  aligned.reference = rep("-",alignment.length)
  aligned.reference[which(exploded.cigar %in% c("M","D","N"))] = strsplit(reference,split = "")[[1]]
  alignment = rbind(aligned.read,aligned.reference)
  colnames(alignment) = exploded.cigar
  alignment.uncut = alignment
  alignment = alignment[,which(colnames(alignment) != "H"),drop = FALSE] # cigar H is for hard-clipping
  alignment = alignment[,which(colnames(alignment) != "S"),drop = FALSE] # cigar S is for soft-clipping
  alignment = alignment[,which(colnames(alignment) != "N"),drop = FALSE] # cigar N is for splicing
  alignment = alignment[,which(alignment["aligned.reference",] != "-"),drop = FALSE] # "-" is for insertions in the reference
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  return(alignment)
}


### update basic objects ###


basic.objects = c(ls(),"alignment.reconstruction")


### alignment.reconstruction ###


bam = get(load(file.path("...","bam.RData")))
aligned.sequence.list.K562 = get(load(file.path("...","aligned.sequence.list.RData")))
cigar.list.K562 = get(load(file.path("...","cigar.list.RData")))
read.to.fast5.name.K562 = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list.K562 = get(load(file.path("...","read.sequence.list.RData")))

registerDoParallel(cores = mc.cores)
alignment.reconstruction = foreach(n = names(bam.K562),.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.sequence.list","read.to.fast5.name"))) %dopar% alignment.reconstruction.function(n)
names(alignment.reconstruction) = names(bam)

save(alignment.reconstruction,file=file.path("...","alignment.reconstruction.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



