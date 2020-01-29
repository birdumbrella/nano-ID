

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


### extract fastq files ###


dir.create(file.path(file.path(prewd,"...")))

for (run in c("...","...")){
  fast5.files = list.files(file.path(prewd,"..."),pattern = ".fast5$",full.names = TRUE,recursive = TRUE)
  
  sink(file.path(prewd,"...",paste0(run,".fastq")))
  for (fast5.file in fast5.files){cat(rna.to.dna(h5read(fast5.file,"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"),suffix = "..."))}
  sink()
}

rm(list = setdiff(ls(),basic.objects))
gc()



