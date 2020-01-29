

### prewd ###


prewd = file.path(getwd(),"...")


### libraries ###


library(foreach)
library(doParallel)
library(rhdf5)


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects","fast5.file","fast5.files")


### load objects ###


load(file.path("...","fast5.files.RData"))

objects = list.files(file.path(prewd,"BasicObjects"))
for (object in objects){source(file.path(prewd,"BasicObjects",object))}


### build.read.speed.assessment.list ###


build.read.speed.assessment.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    
    read.number = strsplit(strsplit(fast5.file,split = "read_")[[1]][2],split = "_")[[1]][1]
    if (is.na(read.number)){read.number = strsplit(strsplit(fast5.file,split = "read")[[1]][3],split = "_")[[1]][1]}
    
    if (is.na(read.number)){move = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Analyses/Basecall_1D_000/BaseCalled_template/Move"))} else {move = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Move")}
    move.rle = Rle(move)

    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    
    results = c(read.name,
                1/(((cumsum(event.repeats[1:1000])/1:1000)*10)/hertz[read.name])
    )
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",1:1000)


### update basic objects ###


basic.objects = c(ls(),"build.read.speed.assessment.list","col.names")


### read speed assessment ###


dir.create(file.path("...","RawSignal"))

bam = get(load(file.path("...","bam.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))
sequencing.summary = get(load(file.path("...","sequencing.summary.RData")))

hertz = rep(...,length(names(bam)))
names(hertz) = names(bam)

index.list = splitvector(1:length(names(bam)),100)

read.speed.assessment.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  read.speed.assessment.subset = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary","hertz"))) %dopar% build.read.speed.assessment.list(n)
  read.speed.assessment.list[[j]] = Reduce('rbind',read.speed.assessment.subset)
}

read.speed.assessment = Reduce('rbind',read.speed.assessment.list)
colnames(read.speed.assessment) = col.names
rownames(read.speed.assessment) = read.speed.assessment[,"read.name"]
read.speed.assessment = apply(read.speed.assessment[,setdiff(col.names,"read.name")],c(1,2),as.numeric)

save(read.speed.assessment,file=file.path("...","read.speed.assessment.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



