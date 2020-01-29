

### prewd ###


prewd = file.path(getwd(),"...")


### libraries ###


library(rhdf5)


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects","fast5.file","fast5.files")


### create fast5.files vector ###


fast5.files = list.files(file.path(prewd,"..."),pattern = ".fast5$",full.names = TRUE,recursive = TRUE)
fast5.files = substr(fast5.files,nchar(prewd)+2,nchar(fast5.files))

save(fast5.files,file=file.path("...","fast5.files.RData"))


### extract read sequences and names ###


dir.create(file.path("..."))

read.to.fast5.name = c()
read.sequence.list = list()
for (fast5.file in fast5.files){
  try({
    h5 = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq")
    read.name = strsplit(strsplit(h5,split = "\n")[[1]][1],split = " ")[[1]][1]
    read.to.fast5.name[paste0(substring(read.name,2,nchar(read.name)),"...")] = fast5.file
    read.sequence.list[[paste0(substring(read.name,2,nchar(read.name)),"...")]] = gsub("U","T",strsplit(h5,split = "\n")[[1]][2])
  },silent = TRUE)
}

save(read.to.fast5.name,file=file.path("...","read.to.fast5.name.RData"))
save(read.sequence.list,file=file.path("...","read.sequence.list.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()



