\name{MVP.PCA}
\alias{MVP.PCA}
\title{Principal Component Analysis}
\description{
Object: Principal Component Analysis
}
\usage{
MVP.PCA (M, perc=1, pcs.keep=5, memo=NULL)
}
\arguments{
	\item{M}{Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size}
	\item{perc}{percentage of total SNPs selected for PCA}
	\item{pcs.keep}{maximum number of PCs}
	\item{memo}{a marker on output file name}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
