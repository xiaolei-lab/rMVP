# MVP [![](https://img.shields.io/badge/Issues-1%2B-brightgreen.svg)](https://github.com/XiaoleiLiuBio/MVP/issues) [![](https://img.shields.io/badge/Release-v1.0.1-blue.svg)](https://github.com/XiaoleiLiuBio/MVP/commits/master)
## A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study
### Installation
**MVP** is only available on GitHub, and can be installed using **devtools**. Two packages should be installed beforehand, **snpStats** and **rfunctions** (only accepts **RcppEigen** <= "0.3.2.9.0"). **MVP** can be installed with the following R code:
```r
#if "devtools" isn't installed, please "install.packages('devtools')" first.
install_version('RcppEigen', version = "0.3.2.9.0")
devtools::install_github("Bioconductor-mirror/snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("xiaoleiLiubio/MVP")
```
After installed successfully, **MVP** can be loaded with the following R code:
```r
library(MVP)
```
Typing ```?MVP``` could get the details of all parameters.

### 1. Data Preparation

1.1 If you have genotype data in **PLINK** format (bed/bim/fam):<br>
**fileBed**, the name of genotype data in PLINK format<br>
**fileKin** is "TRUE" or "FALSE", if true, a kinship matrix represents relationship among individuals will be calculated<br>
**filePC** is "TRUE" or "FALSE", if true, principal component analysis will be performed<br>
**out**, the name of output file<br>
**priority** is "speed" or "memory", the 'speed' mode is faster but uses more memory while 'memory' is slower but uses less memory<br>
**maxLine** is a number, if **priority = "memory"**, it is the number of markers read into memory<br>
```r
MVP.Data(fileBed="plink",
         filePhe=NULL,
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.plink",         
         #priority="memory",
         #maxLine=10000,
         )
```
1.2 If you have genotype data in **Hapmap** format (bed/bim/fam):<br>
**fileHMP** is a string or a string vector, e.g. fileHMP = "hapmap.txt" or fileHMP = c("chr1.hmp.txt", "chr2.hmp.txt", chr3.hmp.txt)<br>
**filePhe**, name of phenotype file<br>
**sep.hmp**, seperator of hapmap file<br>
**sep.phe**, seperator of phenotype file<br>
**SNP.effect** is "Add" or "Dom"<br>
**fileKin** is "TRUE" or "FALSE", if true, a kinship matrix represents relationship among individuals will be calculated<br>
**filePC** is "TRUE" or "FALSE", if true, principal component analysis will be performed<br>
**out**, the name of output file<br>
**priority** is "speed" or "memory", the 'speed' mode is faster but uses more memory while 'memory' is slower but uses less memory<br>
**maxLine** is a number, if **priority = "memory"**, it is the number of markers read into memory<br>
```r
MVP.Data(fileHMP="hapmap.txt",
         filePhe="Phenotype.txt",
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.hmp",
         #priority="memory",
         #maxLine=10000
         )
```
1.3 If you have genotype data in **Numeric** format (bed/bim/fam):<br>
**fileNum**, the name of genotype data in PLINK format<br>
**filePhe**, name of phenotype file<br>
**fileMap**, name of map file<br>
**sep.num**, seperator of Numeric file<br>
**sep.phe**, seperator of phenotype file<br>
**type.geno**, the type of data in Numeric file, "char", "integer", or "double"<br>
**fileKin** is "TRUE" or "FALSE", if true, a kinship matrix represents relationship among individuals will be calculated<br>
**filePC** is "TRUE" or "FALSE", if true, principal component analysis will be performed<br>
**out** is a string, the name of output file<br>
**priority** is "speed" or "memory", the 'speed' mode is faster but uses more memory while 'memory' is slower but uses less memory<br>
**maxLine** is a number, if **priority = "memory"**, it is the number of markers read into memory<br>
```r
MVP.Data(fileNum="Numeric.txt",
         filePhe="Phenotype.txt",
         fileMap="Map.txt",
         sep.num="\t",
         sep.map="\t", 
         sep.phe="\t",
         type.geno="char",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.num",
         #priority="memory"，
         #maxLine=10000
         )
```
1.4 If you have Kinship matrix data that represents the relationship among individuals<br>
**fileKin**, the name of Kinship matrix data, the dimension of Kinship matrix is n * n (n is sample size), no taxa names included<br>
**type.kin**, the type of data in Kinship matrix file, "char", "integer", or "double"<br>
**sep.kin**, seperator of Kinship matrix data file<br>
```r
MVP.Data(fileKin="mvp.kin.txt", 
         type.kin="double",
         sep.kin="\t"
         )
```
1.5 If you have Principal Components matrix data<br>
**filePC**, the name of Principal Components matrix data, the dimension of Principal Components matrix is n * nPC (n is sample size, nPC is number of first columns of PCs), no taxa names and header row included<br>
**type.pc**, the type of data in Principal Components matrix file, "char", "integer", or "double"<br>
**sep.pc**, seperator of Principal Components matrix data file<br>
```r
MVP.Data(filePC="mvp.pc.txt", 
     type.pc="double",
     sep.pc="\t"
     )
```

### 2. Data Input

2.1 At least you should have three data: genotype, phenotype, and map<br>
**genotype**, genotype data generated by **'MVP.Data'** function<br>
**phenotype**, phenotype data, the first column is taxa name and second column is phenotype value<br>
**map**, SNP map information, the first column is SNP name, the second column is Chromosome ID, the third column is phsical position<br>
```r
genotype <- attach.big.matrix("mvp.geno.desc")
phenotype <- read.table("mvp.phe",head=TRUE)
map <- read.table("mvp.map" , head = TRUE)
```
2.2 You can give MVP the prepared Kinship matrix and Covariates data generated by **'MVP.Data'** function<br>
**Kinship**, Kinship matrix, the dimension of Kinship matrix is n * n (n is sample size), no taxa names included<br>
**Covariates**, Covariates matrix, the dimension of Covariates matrix is n * nCV (n is sample size, nCV is number of covariates, no taxa names and header row included<br>
```r
Kinship <- attach.big.matrix("mvp.kin.desc")
Covariates <- attach.big.matrix("mvp.pc.desc")
```
2.3 If you have prepared Kinship matrix and Covariates data generated by other software packages, see **1.4** and **1.5**<br>

### 3. Start GWAS

Three models are included in MVP package: General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU.
**phe**, phenotype data, details see **2.1**<br>
**geno**, genotype data, details see **2.1**<br>
**map**, map data, details see **2.1**<br>
**K**, Kinship matrix, details see **2.2**<br>
**CV.GLM**, Covariates added in GLM, details see **2.2**<br>
**CV.MLM**, Covariates added in MLM, details see **2.2**<br>
**CV.FarmCPU**, Covariates added in FarmCPU, details see **2.2**<br>
**please attention that if nPC.GLM > 0, no PCs should be added in CV.GLM**<br>
**nPC.GLM**, number of first columns of Principal Components added in GLM<br>
**please attention that if nPC.MLM > 0, no PCs should be added in CV.MLM**<br>
**nPC.MLM**, number of first columns of Principal Components added in MLM<br>
**please attention that if nPC.FarmCPU > 0, no PCs should be added in CV.FarmCPU**<br>
**nPC.FarmCPU**, number of first columns of Principal Components added in FarmCPU<br>
**perc**, percentage of random selected SNPs used for calculating Principal Components<br>
**priority**, **"speed"** or **"memory"**<br>
**ncpus**, number of CPUs used for parallel computation<br>
**vc.method**, method of variance components analysis, two methods are avaiblable, "EMMA" and "GEMMA"
```r
iGWAS <- MVP(
    phe=phenotype,
    geno=genotype,
    map=map,
    #K=K,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    perc=1,
    priority="speed",
    ncpus=10,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
    method=c("GLM", "MLM", "FarmCPU")
)
```


