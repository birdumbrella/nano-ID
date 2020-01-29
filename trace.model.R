

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


### build.trace.list ###


build.trace.list = function(read.name){
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
    trace.move = trace[,move > 0]
    trace.no.move = trace[,move == 0]
    
    results = c(read.name,
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_"),seq(0,1,length.out = 10),paste0)))


### update basic objects ###


basic.objects = c(ls(),"build.trace.list","col.names")


### probability of traces ###


dir.create(file.path("..."))

bam = get(load(file.path("...","bam.RData")))
read.to.fast5.name = get(load(file.path("...","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("...","read.sequence.list.RData")))

index.list = splitvector(1:length(names(bam)),100)

trace.list.list = list()

for (j in 1:length(index.list)){
  registerDoParallel(cores = mc.cores)
  trace.list.subset = foreach(n = names(bam)[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.list(n)
  trace.list.list[[j]] = Reduce('rbind',trace.list.subset)
}

traces = Reduce('rbind',trace.list.list)
colnames(traces) = col.names
rownames(traces) = traces[,"read.name"]
traces = apply(traces[,setdiff(col.names,"read.name")],c(1,2),as.numeric)

save(traces,file=file.path("...","traces.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



