\name{MVP.MLM}
\alias{MVP.MLM}
\title{To perform GWAS with MLM model and get the P value of SNPs}
\description{
To perform GWAS with MLM model and get the P value of SNPs
}
\usage{
MVP.MLM(phe, geno, K=NULL, CV=NULL, REML=NULL, priority="speed", genoName=NULL, cpu=2, vc.method="EMMA",
	maxLine=1000, file.output=TRUE, memo="MVP")
}
\arguments{
	\item{phe}{phenotype, n * 2 matrix}
	\item{geno}{genotype, m * n, m is marker size, n is population size}
	\item{CV}{covariates}
	\item{K}{Kinship, Covariance matrix(n * n) for random effects; must be positive semi-definite}
	\item{REML}{a list that contains ve and vg}
	\item{priority}{speed or memory}
	\item{genoName}{name of genotype file, when priority is memory, it is used for attach bigmatrix under windows platform}
	\item{cpu}{number of cpus used for parallel computation}
	\item{vc.method}{the methods for estimating variance component("emma" or "gemma")}
	\item{maxLine}{when the priority is 'memory', users can change this parameter to limit the memory}
	\item{file.output}{whether to output files or not}
	\item{memo}{a marker on output file name}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
