\name{MVP.Data}
\alias{MVP.Data}
\title{To prepare data for MVP package}
\description{
To prepare data for MVP package
}
\usage{
MVP.Data(fileHMP=NULL, fileBed=NULL, fileNum=NULL, filePhe=NULL, fileMap=NULL, fileKin=TRUE,
	filePC=TRUE, out=NULL, sep.hmp="\t", sep.num="\t", sep.map="\t", sep.phe="\t", sep.kin="\t",
	sep.pc="\t", type.geno="char", type.kin="double", type.pc="double", type.map="integer",
	SNP.effect="Add",SNP.impute="Middle", maxLine=10000, maxRecord=1e9, maxInd=1e9, priority="speed", perc=1, pcs.keep=5)
}
\arguments{
	\item{fileHMP}{Genotype in hapmap format}
	\item{fileBed}{Genotype in PLINK binary format}
	\item{fileNum}{Genotype in numeric format; pure 0, 1, 2 matrix; m * n, m is marker size, n is sample size}
	\item{filePhe}{Phenotype, two columns, the first column is taxa name, the second column is trait}
	\item{fileMap}{SNP map information, three columns}{SNP name, Chr, Pos}
	\item{fileKin}{Kinship, n * n matrix, n is sample size}
	\item{filePC}{the PCA file}
	\item{out}{A marker on output file name}
	\item{sep.hmp}{seperator for files}
	\item{sep.num}{seperator for files}
	\item{sep.map}{seperator for files}
	\item{sep.phe}{seperator for files}
	\item{sep.kin}{seperator for files}
	\item{sep.pc}{seperator for files}
	\item{type.geno}{type parameter in bigmemory, genotype data}
	\item{type.kin}{type parameter in bigmemory, Kinship}
	\item{type.pc}{type parameter in bigmemory, PC}
	\item{type.map}{type parameter in bigmemory, map}
	\item{SNP.effect}{"Add" or "Dom"}
	\item{SNP.impute}{"Left", "Middle", "Right"}
	\item{maxLine}{number of SNPs, only used for saving memory when calculate kinship matrix}
	\item{maxRecord}{maximum number for markers}
	\item{maxInd}{maximum number for individuals}
	\item{priority}{"speed" or "memory"}
	\item{perc}{percentage of total SNPs selected for PCA}
	\item{pcs.keep}{maximum number of PCs}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
