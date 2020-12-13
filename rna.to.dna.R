

### rna.to.dna function ###


rna.to.dna = function(h5,suffix = ""){
  fastq = strsplit(h5,split = "\n")
  read.name = strsplit(fastq[[1]][1],split = " ")
  read.name[[1]][1] = paste0(read.name[[1]][1],suffix)
  fastq[[1]][1] = paste(c(unlist(read.name),""),collapse = " ")
  fastq[[1]][2] = gsub("U","T",fastq[[1]][2])
  paste(c(unlist(fastq),""),collapse = "\n")
}



