## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(rMVP)

## ----installation, eval=FALSE----------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("rMVP")

## ----bfile-----------------------------------------------------------------
setwd(tempdir())
out <- "vcf"
vcfPath <- system.file("extdata", "01_vcf", "mvp.vcf", package = "rMVP")
cat(readLines(vcfPath, 15), sep = "\n")
MVP.Data.VCF2MVP(vcfPath, out)
dir()
geno <- attach.big.matrix(paste0(out, ".geno.desc"))
genoInd <- read.table(paste0(out, ".geno.ind"))
map <- read.table(paste0(out, ".map"))

## --------------------------------------------------------------------------
sessionInfo()

