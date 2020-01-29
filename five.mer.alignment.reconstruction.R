

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


### five.mer.alignment.reconstruction.function ###


five.mer.alignment.reconstruction.function = function(read.name){
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
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  aligned.read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.read",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  aligned.reference.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.reference",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  
  exploded.cigar.sequence.five.mers = as.character(cigar.five.mers[as.numeric(stats::filter(as.numeric(c("M" = 1,"I" = 2,"D" = 3,"H" = 4,"S" = 5,"N" = 6)[exploded.cigar]),6^(0:8)[1:5],sides = 1)-(sum(6^(0:8)[1:5])-1))[-c(1:4)]])
  
  five.mers.alignment = rbind(exploded.cigar.sequence.five.mers,
                              "aligned.read.sequence.five.mers" = NA,
                              "aligned.reference.sequence.five.mers" = NA)
  five.mers.alignment["aligned.read.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.read.sequence.five.mers
  five.mers.alignment["aligned.reference.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.reference.sequence.five.mers
  five.mers.alignment = five.mers.alignment[,!(exploded.cigar.sequence.five.mers %in% cigar.five.mers.H.S.N.I.centered)]
  
  return(five.mers.alignment)
}


### update basic objects ###


basic.objects = c(ls(),"five.mer.alignment.reconstruction")


### five.mer.alignment.reconstruction ###


bam = get(load(file.path("...","bam.RData")))
aligned.sequence.list = get(load(file.path("...","aligned.sequence.list.RData")))
cigar.list = get(load(file.path("...","cigar.list.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))

index.list = splitvector(1:length(names(bam)),100)

five.mer.alignment.reconstruction.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  five.mer.alignment.reconstruction.list[[j]] = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.sequence.list","read.to.fast5.name"))) %dopar% five.mer.alignment.reconstruction.function(n)
}

five.mer.alignment.reconstruction = Reduce('c',five.mer.alignment.reconstruction.list)
names(five.mer.alignment.reconstruction) = names(bam)[Reduce('c',index.list)]

save(five.mer.alignment.reconstruction,file=file.path("...","five.mer.alignment.reconstruction.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



