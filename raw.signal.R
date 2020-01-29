

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


### build.raw.signal.list ###


build.raw.signal.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    
    read.number = strsplit(strsplit(fast5.file,split = "read_")[[1]][2],split = "_")[[1]][1]
    if (is.na(read.number)){read.number = strsplit(strsplit(fast5.file,split = "read")[[1]][3],split = "_")[[1]][1]}
    
    if (is.na(read.number)){raw.signal = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Raw/Signal"))} else {raw.signal = h5read(file.path(prewd,fast5.file),paste0("/Raw/Reads/Read_",read.number,"/Signal"))}
    
    if (is.na(read.number)){move = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Analyses/Basecall_1D_000/BaseCalled_template/Move"))} else {move = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Move")}
    move.rle = Rle(move)
    
    read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
    
    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    
    stride = 10
    
    num_events_template = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events_template"]
    num_events = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events"]
    
    anchor = (num_events - num_events_template)*stride
    anchored.raw.signal = raw.signal[anchor:length(raw.signal)]
    
    read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[read.sequence]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
    read.sequence.five.mers = c(NA,NA,read.sequence.five.mers,NA,NA)
    
    events.to.raw = cut(1:length(anchored.raw.signal),(0:num_events_template)*stride,right = FALSE,labels = FALSE)
    event.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = events.to.raw,identity),function(x){anchored.raw.signal[x]})
    raw.means = sapply(event.raw.signal,mean)
    
    cumsum.event.repeats = cumsum(event.repeats)
    
    reevent = cbind(c(1,cumsum.event.repeats[-length(cumsum.event.repeats)] + 1),cumsum.event.repeats)
    colnames(reevent) = c("start","end")

    reevents.to.raw = cut(1:length(anchored.raw.signal),c(((0:num_events_template)*stride)[reevent[,"start"]],num_events_template*stride+1),right = FALSE,labels = FALSE)
    reevent.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = reevents.to.raw,identity),function(x){anchored.raw.signal[x]})
    reevent.raw.means = sapply(reevent.raw.signal,mean)
    
    events = cbind("move" = move,
                   "read.sequence" = rep(read.sequence,times = event.repeats),
                   "read.sequence.five.mers" = rep(read.sequence.five.mers,times = event.repeats))
    
    means = cbind("raw.means" = raw.means,
                  "reevent.raw.means" = rep(reevent.raw.means,times = event.repeats))
    
    events.move = events[move == 1,]
    means.move = means[move == 1,]
    
    events.no.move = events[move == 0,]
    means.no.move = means[move == 0,]
    
    results = c(read.name,
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered,"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.containing,"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.non.T.containing,"reevent.raw.means",drop = FALSE],2,mean),
                
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% setdiff(five.mers.T.containing,five.mers.T.centered),"reevent.raw.means",drop = FALSE],2,mean),
                
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.leading,"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.lagging,"reevent.raw.means",drop = FALSE],2,mean),
                
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered.leading,"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered.lagging,"reevent.raw.means",drop = FALSE],2,mean),
                
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.A.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.C.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.G.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.A.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.C.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.G.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                
                sapply(five.mers.three.mers.centered.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)}),

                sapply(five.mers.three.mers.centered.T.containing.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)}),
                sapply(five.mers.three.mers.centered.non.T.containing.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              "five.mers.T.centered","five.mers.T.containing","five.mers.non.T.containing",
              "five.mers.non.T.centered",
              "five.mers.T.leading","five.mers.T.lagging","five.mers.T.centered.leading","five.mers.T.centered.lagging",
              "five.mers.A.centered.T.containing","five.mers.C.centered.T.containing","five.mers.G.centered.T.containing",
              "five.mers.A.centered.non.T.containing","five.mers.C.centered.non.T.containing","five.mers.G.centered.non.T.containing",
              mkAllStrings(c("A","C","G","T"),3),
              paste0(names(five.mers.three.mers.centered.T.containing.list),".T.containing"),
              paste0(names(five.mers.three.mers.centered.non.T.containing.list),".non.T.containing"))

z.score.rescaling = function(x,y,z){((z - mean(y,na.rm = TRUE))*(sd(x,na.rm = TRUE)/sd(y,na.rm = TRUE))) + mean(x,na.rm = TRUE)}


### update basic objects ###


basic.objects = c(ls(),"build.raw.signal.list","col.names","z.score.rescaling","five.mer.RNA.pore.model.raw.signal","non.T.containing.raw.signal")


### raw signal of kmers ###


dir.create(file.path("..."))

bam = get(load(file.path("...","bam.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))
sequencing.summary = get(load(file.path("...","sequencing.summary.RData")))

index.list = splitvector(1:length(names(bam)),100)

raw.signal.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  raw.signal.subset = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% build.raw.signal.list(n)
  raw.signal.list[[j]] = Reduce('rbind',raw.signal.subset)
}

raw.signal = Reduce('rbind',raw.signal.list)
colnames(raw.signal) = col.names
rownames(raw.signal) = raw.signal[,"read.name"]
raw.signal = apply(raw.signal[,setdiff(col.names,"read.name")],c(1,2),as.numeric)
raw.signal = t(apply(raw.signal,1,function(x){z.score.rescaling(five.mer.RNA.pore.model.raw.signal[names(x[!is.na(x) & non.T.containing.raw.signal])],x[!is.na(x) & non.T.containing.raw.signal],x)}))

save(raw.signal,file=file.path("...","raw.signal.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



