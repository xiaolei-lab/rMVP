\name{MVP-package}
\alias{MVP-package}
\alias{MVP}
\docType{package}
\title{
a Memory efficient, Visualization enhanced, and Parallel accelerated Tool For GWAS
}
\description{
 Object 1: To perform GWAS using General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU model
 Object 2: To calculate kinship among individuals using Varaden method
 Object 3: Estimate variance components using EMMA, FaST-LMM, and HE regression
 Object 4: Generate high-quality figures
}
\usage{
MVP(phe, geno, map, K=NULL, nPC.GLM=NULL, nPC.MLM=NULL, nPC.FarmCPU=NULL,
	perc=1, CV.GLM=NULL,CV.MLM=NULL, CV.FarmCPU=NULL, REML=NULL,
	priority="speed", ncpus=2, vc.method="EMMA",method="MLM", maxLine=1000,
	memo="MVP", P=NULL, method.sub="reward", method.sub.final="reward",
	method.bin="static",bin.size=c(5e5,5e6,5e7),bin.selection=seq(10,100,10),
	Prior=NULL,maxLoop=10, threshold.output=1, iteration.output=FALSE,
	p.threshold=NA, QTN.threshold=NULL,bound=NULL, outward=TRUE,
	col=c("dodgerblue4","olivedrab4","violetred","darkgoldenrod1","purple4"),
	plot.type="b", file.output=TRUE, file="jpg", dpi=300, threshold=0.05,
	Ncluster=1, signal.cex=0.8)
}
\arguments{
	\item{phe}{phenotype, n * 2 matrix, n is sample size}
	\item{geno}{Genotype in bigmatrix format; m * n, m is marker size, n is sample size}
	\item{map}{SNP map information, SNP name, Chr, Pos}
	\item{K}{Kinship, Covariance matrix(n * n) for random effects, must be positive semi-definite}
	\item{nPC.GLM}{number of PCs added as fixed effects in GLM}
	\item{nPC.MLM}{number of PCs added as fixed effects in MLM}
	\item{nPC.FarmCPU}{number of PCs added as fixed effects in FarmCPU}
	\item{perc}{percentage of total SNPs selected for PCA}
	\item{CV.GLM}{covariates added in GLM}
	\item{CV.MLM}{covariates added in MLM}
	\item{CV.FarmCPU}{covariates added in FarmCPU}
	\item{REML}{a list contains ve and vg}
	\item{priority}{speed or memory}
	\item{ncpus}{number of cpus used for parallel}
	\item{vc.method}{methods for estimating variance component("EMMA" or "GEMMA")}
	\item{method}{the GWAS model, "GLM", "MLM", and "FarmCPU", models can be selected simutaneously, i.e. c("GLM", "MLM", "FarmCPU")}
	\item{maxLine}{when the priority is 'memory', users can change this parameter to limit the memory}
	\item{memo}{a marker added on output file name}
	\item{P}{a start p value for each SNP}
	\item{method.sub}{method used in substitution process}
	\item{method.sub.final}{method used in substitution process}
	\item{method.bin}{EMMA or FaSTLMM}
	\item{bin.size}{window size in genome}
	\item{bin.selection}{a vector, how many windows selected}
	\item{Prior}{four columns, SNP name, Chr, Pos, P}
	\item{maxLoop}{maximum number of iterations}
	\item{threshold.output}{output GWAS results only for SNPs with p value lower than the threshold.output}
	\item{iteration.output}{whether to output results for FarmCPU iterations}
	\item{p.threshold}{if all p values in the 1st iteration are bigger than p.threshold, FarmCPU stops}
	\item{QTN.threshold}{Only SNPs have a more significant p value than QTN.threshold have chance to be selected as pseudo QTNs}
	\item{bound}{maximum number of SNPs selected as pseudo QTNs for each iteration}
	\item{outward}{the direction of circular Manhattan plot}
	\item{col}{for color of points in each chromosome on manhattan plot}
	\item{permutation.threshold}{if use a permutation cutoff or not (bonferroni cutoff)}
	\item{permutation.rep}{number of permutation replicates}
	\item{plot.type}{"b" (both Manhattan plot and qq plot will be draw) or "q" (qq plot only)}
	\item{file.output}{whether to output files or not}
	\item{file}{figure format, "jpg", "tiff"}
	\item{dpi}{resolution}
	\item{threshold}{a cutoff line on manhattan plot, 0.05/marker size}
	\item{Ncluster}{number of colors used for drawing PC 1 and PC 2}
	\item{signal.cex}{point size on output figures}
}
\details{
\tabular{ll}{
Package: \tab MVP\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2017-12-12\cr
License: \tab GPL(>=2)\cr
}
}
\author{
Lilin Yin, Haohao Zhang, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
