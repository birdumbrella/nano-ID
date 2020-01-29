

### prewd ###


prewd = file.path(getwd(),"...")


### libraries ###


library(foreach)
library(doParallel)
library(rhdf5)
library(keras)
library(dplyr)


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


### gather data for neural network ###


dir.create(file.path("..."))

load(file.path("...","read.based.mismatch.identification.RData"))
load(file.path("...","traces.RData"))
load(file.path("...","traces.add.on.RData"))
traces = cbind(traces,traces.add.on[rownames(traces),])
load(file.path("...","raw.signal.RData"))
load(file.path("...","raw.signal.five.mers.RData"))
load(file.path("...","raw.signal.five.mers.add.on.RData"))
raw.signal.five.mers = cbind(raw.signal.five.mers,raw.signal.five.mers.add.on[rownames(raw.signal.five.mers),])
read.based.parameters = cbind(read.based.mismatch.identification,traces[rownames(read.based.mismatch.identification),],raw.signal[rownames(read.based.mismatch.identification),],raw.signal.five.mers[rownames(read.based.mismatch.identification),])
colnames(read.based.parameters) = paste0("P",substr(rep("0000",ncol(read.based.parameters)),1,4-nchar(1:ncol(read.based.parameters))),1:ncol(read.based.parameters),"-",colnames(read.based.parameters))

save(read.based.parameters,file=file.path("...","read.based.parameters.RData"))

rm(list = setdiff(ls(),basic.objects))
gc()
  

### neural network ###


load(file.path("...","read.based.parameters.RData"))
read.based.parameters[is.na(read.based.parameters)] = 0

train.model = load_model_hdf5(filepath = file.path("...","NeuralNetwork","train.model.h5"))

dir.create(file.path("..."))

# classification

pred = predict_classes(train.model,read.based.parameters)
names(pred) = rownames(read.based.parameters)

rm(list = setdiff(ls(),basic.objects))
gc()



