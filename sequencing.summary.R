

### prewd ###


prewd = file.path(getwd(),"...")


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects","fast5.file","fast5.files")


### create sequencing summary matrices ###


sequencing.summary = read.delim(file.path(prewd,"...","sequencing_summary.txt"),sep="\t",row.names = NULL,header = TRUE,stringsAsFactors = FALSE)
rownames(sequencing.summary) = sequencing.summary[,"read_id"]

save(sequencing.summary,file=file.path("...","sequencing.summary.RData"))



