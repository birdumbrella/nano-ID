

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


### build.trace.add.on.list ###


build.trace.add.on.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
    trace = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Trace")
    rownames(trace) = c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T")
    move = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Move")
    move.rle = Rle(move)
    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    colnames(trace) = rep(read.sequence,times = event.repeats)
    
    move.rle.A = Rle(move[colnames(trace) == "A"])
    move.rle.C = Rle(move[colnames(trace) == "C"])
    move.rle.G = Rle(move[colnames(trace) == "G"])
    move.rle.T = Rle(move[colnames(trace) == "T"])
    
    results = c(read.name,
                
                mean(runLength(move.rle)[runValue(move.rle) == 1]),
                quantile(runLength(move.rle)[runValue(move.rle) == 1],seq(0,1,length.out = 100)),
                
                mean(runLength(move.rle)[runValue(move.rle) == 0]),
                quantile(runLength(move.rle)[runValue(move.rle) == 0],seq(0,1,length.out = 100)),
                
                tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 0],seq(0,1,length.out = 100))},error = function(x){0})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              paste0("1s","_mean"),
              paste0("1s_",1:100),
              paste0("0s","_mean"),
              paste0("0s_",1:100),
              paste0("1s","_A_mean"),
              paste0("1s_A_",1:100),
              paste0("0s","_A_mean"),
              paste0("0s_A_",1:100),
              paste0("1s","_C_mean"),
              paste0("1s_C_",1:100),
              paste0("0s","_C_mean"),
              paste0("0s_C_",1:100),
              paste0("1s","_G_mean"),
              paste0("1s_G_",1:100),
              paste0("0s","_G_mean"),
              paste0("0s_G_",1:100),
              paste0("1s","_T_mean"),
              paste0("1s_T_",1:100),
              paste0("0s","_T_mean"),
              paste0("0s_T_",1:100))


### update basic objects ###


basic.objects = c(ls(),"build.trace.add.on.list","col.names")


### probability of traces ###


dir.create(file.path("..."))

bam = get(load(file.path("...","bam.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))

index.list = splitvector(1:length(names(bam)),100)

trace.add.on.list.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  trace.add.on.list.subset = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.add.on.list(n)
  trace.add.on.list.list[[j]] = Reduce('rbind',trace.add.on.list.subset)
}

traces.add.on = Reduce('rbind',trace.add.on.list.list)
colnames(traces.add.on) = col.names
rownames(traces.add.on) = traces.add.on[,"read.name"]
traces.add.on = apply(traces.add.on[,setdiff(col.names,"read.name")],c(1,2),as.numeric)

save(traces.add.on,file=file.path("...","traces.add.on.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



