\name{MVP.K.VanRaden}
\alias{MVP.K.VanRaden}
\title{Calculate Kinship matrix by VanRaden method}
\description{
Calculate Kinship matrix by VanRaden method
}
\usage{
MVP.K.VanRaden(M, weight=NULL, priority=c("speed", "memory"), memo=NULL, SUM=NULL, maxLine=1000)
}
\arguments{
	\item{M}{Genotype, m * n, m is marker size, n is population size}
	\item{weight}{vector, the weights for makers}
	\item{priority}{speed or memory}
	\item{memo}{a marker on output file name}
	\item{SUM}{the scaled parameter for Kinship}
	\item{maxLine}{when the priority is 'memory', users can change this parameter to limit the memory}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
