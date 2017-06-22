\name{MVP.FarmCPU}
\alias{MVP.FarmCPU}
\title{Perform GWAS using FarmCPU method}
\description{
Perform GWAS using FarmCPU method
}
\usage{
MVP.FarmCPU(phe, geno, map, CV=NULL, priority="speed", P=NULL, method.sub="reward", method.sub.final="reward",
	method.bin="EMMA", bin.size=c(5e5,5e6,5e7), bin.selection=seq(10,100,10), memo="MVP.FarmCPU",
	Prior=NULL, ncpus=2, maxLoop=10, threshold.output=.01, converge=1, iteration.output=FALSE,
	p.threshold=NA, QTN.threshold=NULL, bound=NULL)
}
\arguments{
	\item{phe}{phenotype, n * 2 matrix}
    \item{phe}{ phenotype}
    \item{geno}{ genotype, n by m matrix. This is Genotype Data Pure (GD). THERE IS NOT COLUMN FOR TAXA.}
    \item{map}{ SNP map information}
    \item{CV}{ covariates}
    \item{priority}{ speed or memory}
    \item{P}{ a start p value for each SNP}
    \item{method.sub, method.sub.final}{ method used in substitution process}
    \item{method.bin}{ EMMA or FaSTLMM}
    \item{bin.size}{ window size in genome}
    \item{bin.selection}{ a vector, how many windows selected}
    \item{memo}{ a marker on output file name}
    \item{Prior}{ four columns, SNP name, Chr, Pos, P}
    \item{ncpus}{ number of cpus used for parallele}
    \item{maxLoop}{ maximum number of iterations}
    \item{threshold.output}{ output GWAS results only for SNPs with p value lower than the threshold.output}
    \item{converge}{ an indicator for when loops stop}
    \item{iteration.output}{ whether to output results for FarmCPU iterations}
    \item{p.threshold}{ if all p values in the 1st iteration are bigger than p.threshold, FarmCPU stops}
    \item{QTN.threshold}{ Only SNPs have a more significant p value than QTN.threshold have chance to be selected as pseudo QTNs}
    \item{bound}{ maximum number of SNPs selected as pseudo QTNs for each iteration}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
