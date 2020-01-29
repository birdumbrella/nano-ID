

### prewd ###


prewd = file.path(getwd(),"...")


### setwd ###


setwd(file.path(prewd,"..."))


### number of cores ###


mc.cores = detectCores()


### basic objects ###


basic.objects = c(ls(),"basic.objects")


### calculate stability based on single molecule counts ###


alpha = log(2)/(24*60)

decay.rate.single.molecule = - alpha - (1/60)*log(1 - labeled.counts/total.counts)
decay.rate.single.molecule[decay.rate.single.molecule <= 0] = NA
synthesis.rate.single.molecule = total.counts*(alpha + decay.rate.single.molecule)
half.lives.single.molecule = log(2)/decay.rate.single.molecule

save(decay.rate.single.molecule,file=file.path("...","decay.rate.single.molecule.RData"))
save(synthesis.rate.single.molecule,file=file.path("...","synthesis.rate.single.molecule.RData"))
save(half.lives.single.molecule,file=file.path("...","half.lives.single.molecule.RData"))



