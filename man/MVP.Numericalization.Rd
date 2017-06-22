\name{MVP.Numericalization}
\alias{MVP.Numericalization}
\title{To convert character SNP genotpe to numerical}
\description{
To convert character SNP genotpe to numerical
}
\usage{
MVP.Numericalization(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE)
}
\arguments{
	\item{x}{genotype in character}
	\item{bit}{number of bit of each genotype}
	\item{effect}{"add", "Dom"}
	\item{impute}{"Middle", "Major", "Minor"}
	\item{Create.indicator}{}
	\item{Major.allele.zero}{}
}
\author{
Lilin Yin, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu \cr
Maintainer: Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}>
}
