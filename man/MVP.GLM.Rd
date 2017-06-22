\name{MVP.GLM}
\alias{MVP.GLM}
\title{To perform GWAS with GLM model and get the P value of SNPs}
\description{
To perform GWAS with GLM model and get the P value of SNPs
}
\usage{
MVP.GLM(phe, geno, CV=NULL, priority="speed", genoName=NULL, cpu=2, memo="MVP.GLM")
}
\arguments{
	\item{phe}{phenotype, n * 2 matrix}
	\item{geno}{genotype, m * n, m is marker size, n is population size}
	\item{CV}{covariates}
	\item{priority}{speed or memory}
	\item{genoName}{name of genotype file, when priority is memory, it is used for attach bigmatrix under windows platform}
	\item{cpu}{number of cpus used for parallel computation}
	\item{memo}{a marker on output file name}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
