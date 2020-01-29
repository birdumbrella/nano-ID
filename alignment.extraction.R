

### prewd ###


prewd = file.path(getwd(),"...")


### libraries ###


library(Biostrings)


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects","fast5.file","fast5.files")


### load objects ###


load(file.path("...","fast5.files.RData"))


### extract aligned sequences ###


reference.genome = readDNAStringSet(filepath = file.path(prewd,"...",".fa"))

dir.create(file.path("..."))

bam = readGAlignments(file = file.path(prewd,"...",".bam"),param = NULL,use.names = TRUE)

aligned.sequence.list = as.list(reference.genome[as(bam,"GRanges")])
aligned.sequence.list = lapply(aligned.sequence.list,as.character)
names(aligned.sequence.list) = names(bam)

cigar.strings = cigar(bam)
cigar.ops = explodeCigarOps(cigar.strings)
cigar.list = explodeCigarOpLengths(cigar.strings)
cigar.list = lapply(1:length(cigar.strings),function(x){vec = cigar.list[[x]];names(vec) = cigar.ops[[x]];vec})
names(cigar.list) = names(bam)

save(bam,file=file.path("...","bam.RData"))
save(aligned.sequence.list,file=file.path("...","aligned.sequence.list.RData"))
save(cigar.list.K562.5EU.0.unlabeled.I,file=file.path("...","cigar.list.RData"))
  
rm(list = setdiff(ls(),basic.objects))
gc()



