# MVP [![](https://img.shields.io/badge/Issues-1%2B-brightgreen.svg)](https://github.com/XiaoleiLiuBio/MVP/issues) [![](https://img.shields.io/badge/Release-v1.0.1-blue.svg)](https://github.com/XiaoleiLiuBio/MVP/commits/master)

## A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/mvp_logo.png">
<img src="results/mvp_logo.png" height="250px" width="450px">
</a>
</p>

### Authors:

> Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, ***Xiaolei Liu***

### Contact:
> xiaoleiliu@mail.hzau.edu.cn

### Contents
* [Installation](#installation)
* [Data Preparation](#data-preparation)
  - [PLINK binary](#plink-binary)
  - [Hapmap](#hapmap)
  - [Numeric](#numeric)
  - [Kinship](#kinship)
  - [Principal Components](#principal-components)
* [Data Input](#data-input)
  - [Basic](#basic)
  - [Advanced](#advanced)
* [Start GWAS](#start-gwas)
* [Output](#output)
  - [SNP-density plot](#SNP-density-plot)
  - [Circular-Manhattan plot](#circular-manhattan-plot)
  - [Rectangular-Manhattan plot](#rectangular-manhattan-plot)
  - [Q-Q plot](#q-q-plot)
* [FAQ AND HINTS](#faq-and-hints)

---
### Installation
**WE STRONGLY RECOMMEND YOU INSTALL MVP ON Microsoft R Open (https://mran.microsoft.com/download/)**<br>
**MVP** is only available on GitHub, and can be installed using **devtools**. Two packages should be installed beforehand, **snpStats** and **rfunctions** (only accepts **RcppEigen** <= "0.3.2.9.0"). **MVP** can be installed with the following R code:<br>
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

---

### Data Preparation

#### PLINK binary
If you have genotype data in **PLINK** format (bed/bim/fam):<br>
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
#### Hapmap
If you have genotype data in **Hapmap** format (bed/bim/fam):<br>
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

| rs# | alleles | chrom | pos | strand | assembly# | center | protLSID | assayLSID | panelLSID | QCcode | A048005080 | A048006063 | A048006555 | A048007096 | A048010273 | ... | A084292044 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| rs3683945 | G/A | 1 | 3197400 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs3707673 | A/G | 1 | 3407393 | + | NA | NA | NA | NA | NA | NA | GA | GA | AA | GA | AA | ... | GG |
| rs6269442 | G/A | 1 | 3492195 | + | NA | NA | NA | NA | NA | NA | AG | GG | GG | AG | GG | ... | AA |
| rs6336442 | G/A | 1 | 3580634 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs13475699 | G | 1 | 3860406 | + | NA | NA | NA | NA | NA | NA | GG | GG | GG | GG | GG | ... | GG |


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
#### Numeric
If you have genotype data in **Numeric** format (bed/bim/fam):<br>
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

| 1 | 1 | 2 | 1 | 2 | … | 0 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: |
| 1 | 1 | 0 | 1 | 0 | … | 2 |
| 1 | 2 | 2 | 1 | 2 | … | 0 |
| 1 | 1 | 2 | 1 | 2 | … | 0 |
| 0 | 0 | 0 | 0 | 0 | … | 0 |

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

#### Kinship
If you have Kinship matrix data that represents the relationship among individuals<br>
**fileKin**, the name of Kinship matrix data, the dimension of Kinship matrix is n * n (n is sample size), no taxa names included<br>
**type.kin**, the type of data in Kinship matrix file, "char", "integer", or "double"<br>
**sep.kin**, seperator of Kinship matrix data file<br>

| 0.3032 | -0.0193 | 0.0094 | 0.0024 | 0.0381 | ... | -0.0072 |
| :---: | :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
| -0.0193 | 0.274 | -0.0243 | 0.0032 | -0.0081 | ... | 0.0056 |
| 0.0094 | -0.0243 | 0.3207 | -0.0071 | -0.0045 | ... | -0.0407 |
| 0.0024 | 0.0032 | -0.0071 | 0.321 | -0.008 | ... | -0.0093 |
| 0.0381 | -0.0081 | -0.0045 | -0.008 | 0.3498 | ... | -0.0238 |
| ... | ... | ... | ... | ... | ... | ... | 
| -0.0072 | 0.0056 | -0.0407 | -0.0093 | -0.0238 | ... | 0.3436 |

```r
MVP.Data(fileKin="mvp.kin.txt", 
         type.kin="double",
         sep.kin="\t"
         )
```

#### Principal Components
If you have Principal Components matrix data<br>
**filePC**, the name of Principal Components matrix data, the dimension of Principal Components matrix is n * nPC (n is sample size, nPC is number of first columns of PCs), no taxa names and header row included<br>
**type.pc**, the type of data in Principal Components matrix file, "char", "integer", or "double"<br>
**sep.pc**, seperator of Principal Components matrix data file<br>
```r
MVP.Data(filePC="mvp.pc.txt", 
     type.pc="double",
     sep.pc="\t"
     )
```

---

### Data Input

#### Basic
At least you should have three data: genotype, phenotype, and map<br>
**genotype**, genotype data generated by **'MVP.Data'** function<br>
**phenotype**, phenotype data, the first column is taxa name and second column is phenotype value<br>
**map**, SNP map information, the first column is SNP name, the second column is Chromosome ID, the third column is phsical position<br>
```r
genotype <- attach.big.matrix("mvp.geno.desc")
phenotype <- read.table("mvp.phe",head=TRUE)
map <- read.table("mvp.map" , head = TRUE)
```

#### Advanced
You can give MVP the prepared Kinship matrix and Covariates data generated by **'MVP.Data'** function<br>
**Kinship**, Kinship matrix, the dimension of Kinship matrix is n * n (n is sample size), no taxa names included<br>
**Covariates**, Covariates matrix, the dimension of Covariates matrix is n * nCV (n is sample size, nCV is number of covariates, no taxa names and header row included<br>
```r
Kinship <- attach.big.matrix("mvp.kin.desc")
Covariates <- attach.big.matrix("mvp.pc.desc")
```
If you have prepared Kinship matrix and Covariates data generated by other software packages, see **Kinship**[#kinship] and **Principal Components**[#principal-components]<br>

---

### Start GWAS

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
**vc.method**, method of variance components analysis, two methods are avaiblable, "EMMA" (Kang, 2008, Genetics) and "GEMMA" (HE Regression, Xiang Zhou, 2016, doi: https://doi.org/10.1101/042846)<br>
**maxLoop**, a parameter for FarmCPU only, the maximum iterations allowed in FarmCPU<br>
**method.bin**, a parameter for FarmCPU only, there are three options: "FaST-LMM","EMMA", and "static"<br>
**method**, models for association tests, three models are available in MVP, **"GLM"**, **"MLM"**, and **"FarmCPU"**, one or two or three models can be selected for assocation tests<br>
```r
imMVP <- MVP(
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

---

### Output
**MVP** will automatically output one(".jpg" by default) of three types (".jpg",".pdf",".tiff") high quality visual plot. Still, users could adjust about 40 parameters to plot more elaborate results by the function of ```MVP.Report()```. 
Generally, ```MVP.Report()``` could accept the final return of ```MVP()``` directly, for example ```MVP.Report(imMVP, ...)```. Nevertheless, users could load the prepared data into R and use the function to visualize it, if in this case, the data at least contains four columns, which are names of SNP, chromosome, postion and P-value of a trait respectively, more traits could be sequentially appended after the data by column. Typing ```?MVP.Report()``` to see the details of all parameters. Typing ```data(pig60K); data(cattle50K)``` to load the attached datasets.

#### SNP-density plot

```r
MVP.Report(pig60K[, c(1:3)], plot.type="d", col=c("darkgreen", "yellow", "red"), file="jpg", dpi=300)
```

**plot.type**, which type will be plotted, if "d", plotting ***SNP-density plot***; if "c", plotting ***Circular-Manhattan plot***; if "m", plotting ***Rectangular-Manhattan plot***, if "q", plotting ***QQ-plot***<br>
**bin.size**, the window size for counting the number of SNP<br>
**bin.max**, the max value of legend, the windows whose number of SNP are bigger than **bin.max** will be painted in same color

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/illumilla_60K.jpg">
<img src="results/illumilla_60K.jpg" height="460px" width="680px">
</a>
</p>

#### Circular-Manhattan plot

***Genome-wide association study(GWAS)***

```r
MVP.Report(pig60K, plot.type="c", chr.labels=paste("Chr",c(1:18,"X"),sep=""), threshold=c(0.05,0.01),
      cir.chr.h=1, amplify=TRUE, threshold.lty=c(2,1), threshold.col=c("blue","red"), signal.line=1,
      signal.col="red", file="jpg", dpi=300)
```

**chr.labels**, renamed names of each chromosome<br>
**threshold**, the significant level for Bonferroni adjustment<br>
**cir.chr.h**, the width of outer circle<br>
**amplify**, "TRUE" or "FALSE", whether to highlight the significant SNPs<br>
**signal.line**, the width of the lines that cross all circle, if signal.line=NULL, the lines that crosse circles won't be added<br>
**signal.col**, the color for the significant SNPs, if NULL, it will use the color index of ***col***<br>
**signal.cex**, the cex for the significant SNPs<br>

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/Circular-Manhattan.jpg">
<img src="results/Circular-Manhattan.jpg" height="400px" width="400px">
</a>
</p>

***Genomic Selection/Prediction(GS/GP)***

```r
MVP.Report(cattle50K, plot.type="c", LOG10=FALSE, outward=TRUE, chr.labels=paste("Chr",c(1:29),sep=""),
         r=1.2, cir.chr.h=1.3, cir.legend.cex=0.5, cir.band=1, threshold=NULL, file="jpg", dpi=300)
```

**LOG10**, "TRUE" or "FALSE", if FALSE, the original values will be used to plot<br>
**outward**, "TRUE" or "FALSE", the direction of plotting points<br>
**r**, the radius of the circle<br>
**cir.legend.cex**, the size of axis number of legend<br>
**cir.band**, the size of interval between circles<br>

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/Circular-Manhattan.cattle.jpg">
<img src="results/Circular-Manhattan.cattle.jpg" height="400px" width="400px">
</a>
</p>

#### Rectangular-Manhattan plot

***Genome-wide association study(GWAS)***

```r
MVP.Report(pig60K[,c(1:3,6)], plot.type="m", threshold=NULL, file="jpg", dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/Rectangular-Manhattan.trait3.jpg">
<img src="results/Rectangular-Manhattan.trait3.jpg" height="300px" width="900px">
</a>
</p>

***Genomic Selection/Prediction(GS/GP)***

```r
MVP.Report(cattle50K[,c(1:3,5)], plot.type="m", LOG10=FALSE, ylab="SNP effect", 
         threshold=NULL, file="jpg", dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/Rectangular-Manhattan.Fat percentage.jpg">
<img src="results/Rectangular-Manhattan.Fat percentage.jpg" height="300px" width="900px">
</a>
</p>

***Multiple tracks rectangular-Manhattan plot***

```r
MVP.Report(imMVP, plot.type="m", threshold=0.05, multracks=TRUE, file="jpg", dpi=300)
```
**multracks**, logical, if FALSE, plotting multiple rectangular Manhattan plots on multiple tracks, if TRUE, all Manhattan plots will be plotted in only one track<br>

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/MVP.Multracks.Rectangular_Manhattan.trait.GLM.trait.MLM.trait.FarmCPU.jpg">
<img src="results/MVP.Multracks.Rectangular_Manhattan.trait.GLM.trait.MLM.trait.FarmCPU.jpg" height="600px" width="500px">
</a>
</p>

#### Q-Q plot

***Single track Q-Q plot***

```r
MVP.Report(pig60K[,c(1:3,6)], plot.type="q", conf.int=TRUE, conf.int.col="grey", file="jpg", dpi=300)
```
**conf.int**, logical, whether to drew the confidence interval on QQ-plot<br>
**conf.int.col**, character, the color of the confidence interval on QQ-plot<br>

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/QQplot.trait3.jpg">
<img src="results/QQplot.trait3.jpg" height="400px" width="400px">
</a>
</p>

***Multiple track Q-Q plot***

```r
MVP.Report(imMVP, plot.type="q", multracks=TRUE, conf.int=TRUE, conf.int.col="grey", file="jpg", dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/XiaoleiLiuBio/MVP/master/results/MVP.Multraits.QQplot.jpg">
<img src="results/MVP.Multraits.QQplot.jpg" height="500px" width="500px">
</a>
</p>

---

### FAQ and Hints

:sos: **Question1:** Failing to install "devtools":

***ERROR: configuration failed for package ‘git2r’***

***removing ‘/Users/acer/R/3.4/library/git2r’***

***ERROR: dependency ‘git2r’ is not available for package ‘devtools’***

***removing ‘/Users/acer/R/3.4/library/devtools’***

:yum: **Answer:** Please type the following codes in terminal.
```ssh
apt-get install libssl-dev/unstable
```
---
:sos: **Question2:** When installing packages from Github with "devtools", there is a error:
 
 ***Error in curl::curl_fetch_disk(url, x$path, handle = handle): Problem with the SSL CA cert (path? access rights?)***
 
:yum: **Answer:** Please type the following codes and than try agian.
```r
library(httr)
set_config(config(ssl_verifypeer = 0L))
```

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/XiaoleiLiuBio/MVP/issues)
