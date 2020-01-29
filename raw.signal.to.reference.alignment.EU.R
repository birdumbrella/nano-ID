

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


### build.signal.position.list ###


build.signal.position.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    
    read.number = strsplit(strsplit(fast5.file,split = "read_")[[1]][2],split = "_")[[1]][1]
    if (is.na(read.number)){read.number = strsplit(strsplit(fast5.file,split = "read")[[1]][3],split = "_")[[1]][1]}
    
    raw.signal = h5read(file.path(prewd,fast5.file),paste0("/Raw/Reads/Read_",read.number,"/Signal"))
    
    move = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_001/BaseCalled_template/Move")
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
                   "read.sequence.five.mers" = rep(read.sequence.five.mers,times = event.repeats),
                   "raw.means" = raw.means,
                   "reevent.raw.means" = rep(reevent.raw.means,times = event.repeats))
    
    read = read.sequence.list[[read.name]]
    reference = aligned.sequence.list[[read.name]]
    cigar = cigar.list[[read.name]]
    exploded.cigar = as.vector(Rle(names(cigar),cigar))
    alignment.length = length(exploded.cigar)
    aligned.read = rep("-",alignment.length)
    aligned.read[which(exploded.cigar %in% c("M","I","S"))] = strsplit(read,split = "")[[1]]
    aligned.reference = rep("-",alignment.length)
    aligned.reference[which(exploded.cigar %in% c("M","D"))] = strsplit(reference,split = "")[[1]]
    alignment = rbind(aligned.read,aligned.reference)
    colnames(alignment) = exploded.cigar
    
    aligned.read.sequence.five.mers = unlist(lapply(1:(nchar(paste(alignment["aligned.read",which(colnames(alignment) != "S")],collapse = "")) - 4),function(i){substr(paste(alignment["aligned.read",which(colnames(alignment) != "S")],collapse = ""),i,5+i-1)}))
    aligned.reference.sequence.five.mers = unlist(lapply(1:(nchar(paste(alignment["aligned.reference",which(colnames(alignment) != "S")],collapse = "")) - 4),function(i){substr(paste(alignment["aligned.reference",which(colnames(alignment) != "S")],collapse = ""),i,5+i-1)}))
    
    spikein.genome = readDNAStringSet(filepath = file.path(prewd,"RawData","SpikeIns","spikeins.corrected.fa"))
    spikein.sequence = as.vector(spikein.genome["chrS8"])
    spikein.sequence.five.mers = unlist(lapply(1:(nchar(spikein.sequence) - 4),function(i){substr(spikein.sequence,i,5+i-1)}))
    names(spikein.sequence.five.mers) = 3:(length(spikein.sequence.five.mers) + 2)
    
    exploded.cigar.sequence = paste(exploded.cigar,collapse = "")
    exploded.cigar.sequence.five.mers = unlist(lapply(1:(nchar(exploded.cigar.sequence) - 4),function(i){substr(exploded.cigar.sequence,i,5+i-1)}))
    
    cigar.five.mers = mkAllStrings(c("M","I","D","S"),5)
    names(cigar.five.mers) = 1:length(cigar.five.mers)
    cigar.five.mer.mat = t(sapply(strsplit(cigar.five.mers,split = ""),c))
    rownames(cigar.five.mer.mat) = cigar.five.mers
    
    cigar.five.mers.I.centered = rownames(cigar.five.mer.mat[apply(cigar.five.mer.mat,1,function(x){x[3] == "I"}),])
    cigar.five.mers.D.centered = rownames(cigar.five.mer.mat[apply(cigar.five.mer.mat,1,function(x){x[3] == "D"}),])
    cigar.five.mers.S.containing = rownames(cigar.five.mer.mat[apply(cigar.five.mer.mat,1,function(x){any(x == "S")}),])
    
    subsel = !is.na(events[,"read.sequence.five.mers"])
    events = events[subsel,]
    
    reevent.raw.signal = reevent.raw.signal[unique(rep(1:length(event.repeats),times = event.repeats)[subsel])]
    reevent.raw.means = reevent.raw.means[unique(rep(1:length(event.repeats),times = event.repeats)[subsel])]
    
    rle.read.sequence.five.mers = Rle(cumsum(events[,"move"]))
    run.length.rle.read.sequence.five.mers = runLength(rle.read.sequence.five.mers)
    run.value.rle.read.sequence.five.mers = events[,"read.sequence.five.mers"][events[,"move"] != 0]
    cumsum.run.length.rle.read.sequence.five.mers = cumsum(run.length.rle.read.sequence.five.mers)
    
    five.mers.alignment = rbind(exploded.cigar.sequence.five.mers,
                                "aligned.read.sequence.five.mers" = NA,
                                "aligned.reference.sequence.five.mers" = NA,
                                "spikein.sequence.five.mers" = NA,
                                "rev(reverse(run.value.rle.read.sequence.five.mers))" = NA,
                                "run.value.rle.read.sequence.five.mers.number" = NA,
                                "spikein.sequence.reference.position" = NA)
    five.mers.alignment["aligned.read.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.read.sequence.five.mers
    five.mers.alignment["aligned.reference.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.reference.sequence.five.mers
    five.mers.alignment["spikein.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing,cigar.five.mers.I.centered))][1:length((start(bam[read.name]) + 2):(end(bam[read.name]) - 2))] = spikein.sequence.five.mers[as.character((start(bam[read.name]) + 2):(end(bam[read.name]) - 2))]
    five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",which(!(exploded.cigar.sequence.five.mers %in% cigar.five.mers.D.centered))[cumsum(c(1,rev(events[,"move"][events[,"move"] != 0])))[-(sum(events[,"move"] != 0) + 1)]]] = rev(sapply(run.value.rle.read.sequence.five.mers,reverse))
    five.mers.alignment["run.value.rle.read.sequence.five.mers.number",which(!(exploded.cigar.sequence.five.mers %in% cigar.five.mers.D.centered))[cumsum(c(1,rev(events[,"move"][events[,"move"] != 0])))[-(sum(events[,"move"] != 0) + 1)]]] = length(run.value.rle.read.sequence.five.mers):1
    five.mers.alignment["spikein.sequence.reference.position",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing,cigar.five.mers.I.centered))][1:length((start(bam[read.name]) + 2):(end(bam[read.name]) - 2))] = (start(bam[read.name]) + 2):(end(bam[read.name]) - 2)
    five.mers.alignment = five.mers.alignment[,which(five.mers.alignment["spikein.sequence.reference.position",] == limits[1]-5):which(five.mers.alignment["spikein.sequence.reference.position",] == limits[2]+5)]
    five.mers.alignment
    dim(five.mers.alignment)
    
    event.numbers = five.mers.alignment["run.value.rle.read.sequence.five.mers.number",which(five.mers.alignment["spikein.sequence.reference.position",] != "NA" & five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",] != "NA")]
    event.positions = as.numeric(five.mers.alignment["spikein.sequence.reference.position",which(five.mers.alignment["spikein.sequence.reference.position",] != "NA" & five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",] != "NA")])
    lengths = sapply(reevent.raw.signal[as.numeric(event.numbers)],length)
    
    events.for.scaling = five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",which(five.mers.alignment["spikein.sequence.reference.position",] != "NA" & five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",] != "NA")] == five.mers.alignment["spikein.sequence.five.mers",which(five.mers.alignment["spikein.sequence.reference.position",] != "NA" & five.mers.alignment["rev(reverse(run.value.rle.read.sequence.five.mers))",] != "NA")]

    if (sum(events.for.scaling) < 10){events.for.scaling = TRUE = break}
    signal.position.vec = z.score.rescaling(five.mer.RNA.pore.model[run.value.rle.read.sequence.five.mers[as.numeric(event.numbers[events.for.scaling])],"mean"],reevent.raw.means[as.numeric(event.numbers[events.for.scaling])],unlist(lapply(reevent.raw.signal,rev)[as.numeric(event.numbers)]))
    names(signal.position.vec) = unlist(sapply(1:length(event.positions),function(x){seq(event.positions[x]-0.5,event.positions[x]+0.5,length.out = lengths[x])}))
    
    return(signal.position.vec)
  },silent = TRUE)
}
  
z.score.rescaling = function(x,y,z){((z - mean(y,na.rm = TRUE))*(sd(x,na.rm = TRUE)/sd(y,na.rm = TRUE))) + mean(x,na.rm = TRUE)}


### update basic objects ###


basic.objects = c(ls(),"build.signal.position.list","z.score.rescaling")


### raw signal to reference alignment ###


load(file.path("...","bam.RData"))
load(file.path("...","aligned.sequence.list.RData"))
load(file.path("...","cigar.list.RData"))
load(file.path("...","read.to.fast5.name.RData"))
load(file.path("...","read.sequence.list.RData"))
load(file.path("...","sequencing.summary.RData"))

spikein.genome = readDNAStringSet(filepath = file.path(prewd,"...",".fa"))
spikein.sequence = as.vector(spikein.genome)
spikein.sequence.five.mers = unlist(lapply(1:(nchar(spikein.sequence) - 4),function(i){substr(spikein.sequence,i,5+i-1)}))
names(spikein.sequence.five.mers) = 3:(length(spikein.sequence.five.mers) + 2)

model.means = five.mer.RNA.pore.model[sapply(spikein.sequence.five.mers,reverse),"mean"]
model.sds = five.mer.RNA.pore.model[sapply(spikein.sequence.five.mers,reverse),"sd"]
model.positions = as.numeric(names(spikein.sequence.five.mers))

limits = ...

read.sequence.length = sapply(read.sequence.list,nchar)
reads = names(read.sequence.list)[order(read.sequence.length,decreasing = TRUE)]
reads.subsel = intersect(reads,names(subset(bam,start < limits[1]-5 & end > limits[2]+5)))
length(reads.subsel)

registerDoParallel(cores = mc.cores)
signal.position.list = foreach(n = reads.subsel,.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% build.signal.position.list(n)
names(signal.position.list) = reads.subsel

par(mfrow = c(2,1))

signal.position.list.subsel = signal.position.list[1:limit]
x = as.numeric(unlist(lapply(signal.position.list.subsel[sapply(signal.position.list.subsel,class) != "try-error"],names)))
y = as.numeric(unlist(signal.position.list.subsel[sapply(signal.position.list.subsel,class) != "try-error"]))

plot.new()
plot.window(xlim = limits,ylim = range(pretty(quantile(five.mer.RNA.pore.model[,"mean"],probs = c(0.001,0.999),na.rm = TRUE))))
axis(1,lwd = 1)
axis(2,lwd = 1)
box()
title(ylab = "Signal [pA]",xlab = "Position [bp]")
for (i in which(spikein.sequence.five.mers %in% five.mers.T.containing)+2){
  polygon(c(i-0.5,i-0.5,i+0.5,i+0.5),c(40,140,140,40),col = convertcolor("darkgrey",50),lty=0)
}
for (i in which(spikein.sequence.five.mers %in% five.mers.T.centered)+2){
  polygon(c(i-0.5,i-0.5,i+0.5,i+0.5),c(40,140,140,40),col = convertcolor("darkgrey",50),lty=0)
}
abline(v = ((limits[1]-1):limits[2])+0.5,lty = 2,col = "lightgrey")
for (i in 2:length(model.positions)){
  polygon(c(model.positions[i-1]-0.5,model.positions[i-1]-0.5,model.positions[i-1]+0.5,model.positions[i-1]+0.5),c((model.means - model.sds)[i],(model.means + model.sds)[i],(model.means + model.sds)[i],(model.means - model.sds)[i]),col = convertcolor("darkblue",50),lty=0)
  lines(c(model.positions[i-1]-0.5,model.positions[i-1]+0.5),c(model.means[i],model.means[i]),col = "darkblue")
}
points(x,y,pch = 19,cex = 0.1,col = convertcolor("black","05"))
text(1:length(strsplit(spikein.sequence,split = "")[[1]]),rep(50,length(strsplit(spikein.sequence,split = "")[[1]])),strsplit(spikein.sequence,split = "")[[1]])

plot.new()
plot.window(xlim = limits,ylim = range(pretty(quantile(five.mer.RNA.pore.model[,"mean"],probs = c(0.001,0.999),na.rm = TRUE))))
axis(1,lwd = 1)
axis(2,lwd = 1)
box()
title(ylab = "Signal [pA]",xlab = "Position [bp]")
for (i in which(spikein.sequence.five.mers %in% five.mers.T.containing)+2){
  polygon(c(i-0.5,i-0.5,i+0.5,i+0.5),c(40,140,140,40),col = convertcolor("darkgrey",50),lty=0)
}
for (i in which(spikein.sequence.five.mers %in% five.mers.T.centered)+2){
  polygon(c(i-0.5,i-0.5,i+0.5,i+0.5),c(40,140,140,40),col = convertcolor("darkgrey",50),lty=0)
}
abline(v = ((limits[1]-1):limits[2])+0.5,lty = 2,col = "lightgrey")
for (i in 2:length(model.positions)){
  polygon(c(model.positions[i-1]-0.5,model.positions[i-1]-0.5,model.positions[i-1]+0.5,model.positions[i-1]+0.5),c((model.means - model.sds)[i],(model.means + model.sds)[i],(model.means + model.sds)[i],(model.means - model.sds)[i]),col = convertcolor("darkblue",50),lty=0)
  lines(c(model.positions[i-1]-0.5,model.positions[i-1]+0.5),c(model.means[i],model.means[i]),col = "darkblue")
}
split.seq = seq(limits[1]-5,limits[2]+5,0.1)
aggregated = sapply(lapply(tapply(1:length(y),INDEX = cut(x,split.seq,right = FALSE,labels = FALSE),identity),function(x){y[x]}),median)
lines(split.seq[-length(split.seq)] + 0.05,aggregated,lwd = 2)
text(1:length(strsplit(spikein.sequence,split = "")[[1]]),rep(50,length(strsplit(spikein.sequence,split = "")[[1]])),strsplit(spikein.sequence,split = "")[[1]])

rm(list = setdiff(ls(),basic.objects))
gc()



