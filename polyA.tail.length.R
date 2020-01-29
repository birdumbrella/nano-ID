

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


### build.polyA.tail.features.list ###


build.polyA.tail.features.list = function(read.name){
  fast5.file = read.to.fast5.name[read.name]
  
  read.number = strsplit(strsplit(fast5.file,split = "read_")[[1]][2],split = "_")[[1]][1]
  if (is.na(read.number)){read.number = strsplit(strsplit(fast5.file,split = "read")[[1]][3],split = "_")[[1]][1]}
  
  if (is.na(read.number)){raw.signal = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Raw/Signal"))} else {raw.signal = h5read(file.path(prewd,fast5.file),paste0("/Raw/Reads/Read_",read.number,"/Signal"))}
  
  results = c(read.name,length(raw.signal),NA,NA,NA,NA,NA,NA)
  
  try({
    upper.limit = quantile(raw.signal,0.999)
    lower.limit = quantile(raw.signal,0.001)
    
    raw.signal[which(raw.signal > upper.limit | raw.signal < lower.limit)] = stats::filter(raw.signal,rep(1/10,10))[which(raw.signal > upper.limit | raw.signal < lower.limit)]
    raw.signal[is.na(raw.signal)] = upper.limit
    
    subsel.length = length(raw.signal)
    raw.signal.subsel = raw.signal[1:subsel.length]
    trend.lines = sort(kmeans(raw.signal.subsel[which(raw.signal.subsel > lower.limit | raw.signal.subsel < upper.limit)],2)$centers)
    
    upper.trend = (trend.lines[1] + (trend.lines[2] - trend.lines[1])*(60/100))
    lower.trend = (trend.lines[1] + (trend.lines[2] - trend.lines[1])*(40/100))
    
    raw.signal = c(rep(trend.lines[2],5000),raw.signal[-(1:2000)])
    raw.signal.subsel = c(rep(trend.lines[2],5000),raw.signal.subsel[-(1:2000)])
    
    raw.signal.squished = pmin(raw.signal.subsel,upper.trend)
    raw.signal.squished = pmax(raw.signal.squished,seq(lower.trend,upper.trend,length.out = length(raw.signal.squished)))
    
    upper.line.loss = (raw.signal.squished - upper.trend)^2
    lower.line.loss = (raw.signal.squished - seq(lower.trend,upper.trend,length.out = length(raw.signal.squished)))^2
    initial.sum = sum(upper.line.loss)
    cumsum.upper.line.loss = cumsum(upper.line.loss)
    cumsum.lower.line.loss = cumsum(lower.line.loss)
    loss.score = sapply(1:length(raw.signal.squished),function(x){initial.sum - cumsum.upper.line.loss[x] + cumsum.lower.line.loss[x]})
    
    polyA.tail.start = which.min(loss.score)
    adapter.segment.start = which.max(loss.score[1:polyA.tail.start])
    
    window = 1500
    
    sign.indicator = Rle(sign(loss.score[1:(subsel.length - window)] - loss.score[(1:(subsel.length - window)) + window]))
    
    if (length(runValue(sign.indicator)) > 2){
      if (runValue(sign.indicator)[1] == -1){
        polyA.tail.start = which.min(loss.score[runLength(sign.indicator)[1]:sum(runLength(sign.indicator)[1:3])]) + runLength(sign.indicator)[1] - 1
      } else {
        polyA.tail.start = which.min(loss.score)
      }
    } else {
      polyA.tail.start = which.min(loss.score)
    }
    adapter.segment.start = which.max(loss.score[1:polyA.tail.start])
    
    raw.signal.clipped = raw.signal[polyA.tail.start:min(polyA.tail.start + (hertz[read.name]*5000/70),length(raw.signal))]
    
    cumulative.mean = cumsum(raw.signal.clipped)/(1:length(raw.signal.clipped))
    cumulative.standard.deviation = sqrt(cumsum((raw.signal.clipped - cumulative.mean)^2)/((1:length(raw.signal.clipped)) - 1))
    cumulative.standard.deviation[1:100] = NA
    cumulative.standard.deviation[1:(which.min(cumulative.standard.deviation)/2)] = NA
    cumulative.standard.deviation = cumulative.standard.deviation[1:which.max(cumulative.standard.deviation > sum(range(cumulative.standard.deviation,na.rm = TRUE))/2)]
    
    cumulative.standard.deviation = cumulative.standard.deviation - min(cumulative.standard.deviation,na.rm = TRUE)
    
    polyA.tail.signal.length = which.min(cumulative.standard.deviation*length(cumulative.standard.deviation)/max(cumulative.standard.deviation,na.rm = TRUE) - 1:length(cumulative.standard.deviation))
    
    read.sequence.length = nchar(read.sequence.list[[read.name]])
    read.signal.length = length(raw.signal) - (polyA.tail.start + polyA.tail.signal.length)
    
    read.speed.bp.sec = read.sequence.length/(read.signal.length/hertz[read.name])
    
    polyA.tail.length = polyA.tail.signal.length/hertz[read.name] + 5
    corrected.polyA.tail.length = read.speed.bp.sec*polyA.tail.signal.length/hertz[read.name] + 5
    
    results = c(read.name,length(raw.signal),read.sequence.length,read.speed.bp.sec,polyA.tail.signal.length,hertz[read.name],polyA.tail.length,corrected.polyA.tail.length)
    names(results) = NULL
  },silent = TRUE)
  
  return(results)
}

col.names = c("read.name","raw.signal.length","read.sequence.length","read.speed.bp.sec","polyA.tail.signal.length","hertz","polyA.tail.length","corrected.polyA.tail.length")


### update basic objects ###


basic.objects = c(ls(),"build.polyA.tail.features.list","col.names")


### polyA-tail length ###


dir.create(file.path("..."))

bam = get(load(file.path("...","bam.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))
sequencing.summary = get(load(file.path("...","sequencing.summary.RData")))

hertz = rep(...,length(names(bam)))
names(hertz) = names(bam)

index.list = splitvector(1:length(names(bam)),100)

polyA.tail.features.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  polyA.tail.features.subset = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary","hertz"))) %dopar% build.polyA.tail.features.list(n)
  polyA.tail.features.list[[j]] = Reduce('rbind',polyA.tail.features.subset)
}

polyA.tail.features = Reduce('rbind',polyA.tail.features.list)
colnames(polyA.tail.features) = col.names
rownames(polyA.tail.features) = polyA.tail.features[,"read.name"]
polyA.tail.features = apply(polyA.tail.features[,setdiff(col.names,"read.name")],c(1,2),as.numeric)

save(polyA.tail.features,file=file.path("...","polyA.tail.features.RData"))
  
rm(list = setdiff(ls(),basic.objects))
gc()



