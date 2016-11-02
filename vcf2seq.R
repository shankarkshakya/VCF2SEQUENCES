rm(list=ls())

library(vcfR) 
library(ape)

all_PITG <- list.files("Bedified_INF_IPO_MIR/") #list of vcfs for each PITG
length(all_PITG)

pitg_cov <- read.csv("pitg_coverage.csv", header = TRUE) # this file has coordinates for each PITG


dna <- ape::read.dna(file = "phytophthora_infestans_t30-4_1_genes.fasta", format="fasta") # this file has PITG sequences and will be used downstream to generate reference sequence.


trees <- vector("list", 20)
nt_div <- vector("list", 20)


for (k in 1:length){
  PITG_vcf <- read.vcfR(file.path("Bedified_INF_IPO_MIR/", all_PITG[k]))
  PITG_vcf
  
  #dna <- ape::read.dna(file = "phytophthora_infestans_t30-4_1_genes.fasta", format="fasta")
  
  seq <- dna[k]
  #seq <- dna[ grep("PITG_00002", names(dna)) ]
  
  my_dnabin <- vcfR2DNAbin(PITG_vcf, consensus = FALSE, extract.haps = TRUE, gt.split = "/", 
                           ref.seq = seq, start.pos = pitg_cov[k,3])
  
  my_dnabin
  
  #checkAlignment(my_dnabin)
  
  #write.dna(my_dnabin, file = paste(all_PITG[k], "fasta", sep = ".") , format = "fasta")
  
  
  #ape::image.DNAbin(my_dnabin[,ape::seg.sites(my_dnabin)])
  
  
  my_dist <- dist.dna(my_dnabin, variance = TRUE)
  
  trees[[k]] <- nj(my_dist)
  
  nt_div[[k]] <- nuc.div(my_dnabin)
  plot(root(nj(my_dist), outgroup = "mir1_0"), type = "phylogram")
  #add.scale.bar(cex = 0.8, font = 2, col = "black")
  
  
}


```


## Looping through list of pitgs to convert from vcf to dnabin

```{r}

pitglist #list of pitg you want to convert from vcf to dna bin

trees <- vector("list", length(pitglist))
nt_div <- vector("list", length(pitglist))


for (j in 1:length(pitglist)){
  k <- grep(pitglist[j], all_PITG)
  
  PITG_vcf <- read.vcfR(file.path("Bedified_INF_IPO_MIR/", all_PITG[k]))
  PITG_vcf
  
  #dna <- ape::read.dna(file = "phytophthora_infestans_t30-4_1_genes.fasta", format="fasta")
  
  seq <- dna[k]
  #seq <- dna[ grep("PITG_00002", names(dna)) ]
  
  my_dnabin <- vcfR2DNAbin(PITG_vcf, consensus = FALSE, extract.haps = TRUE, gt.split = "/", 
                           ref.seq = seq, start.pos = pitg_cov[k,3])
  
  my_dnabin
  
  #checkAlignment(my_dnabin)
  
  #write.dna(my_dnabin, file = paste(all_PITG[k], "fasta", sep = ".") , format = "fasta")
  
  
  #ape::image.DNAbin(my_dnabin[,ape::seg.sites(my_dnabin)])
  
  
  my_dist <- dist.dna(my_dnabin, variance = TRUE)
  
  trees[[j]] <- nj(my_dist)
  
  nt_div[[j]] <- nuc.div(my_dnabin)
}

names(trees) <- pitglist
