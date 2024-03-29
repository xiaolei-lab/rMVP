% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MVP.BRENT.Vg.Ve.R
\name{MVP.BRENT.Vg.Ve}
\alias{MVP.BRENT.Vg.Ve}
\title{MVP.BRENT.Vg.Ve
variance component estimation using the BRENT method}
\usage{
MVP.BRENT.Vg.Ve(y, X, eigenK, verbose = FALSE)
}
\arguments{
\item{y}{phenotype}

\item{X}{covariate matrix, the first column is 1s}

\item{eigenK}{eigen of Kinship matrix}

\item{verbose}{whether to print detail.}
}
\value{
vg, ve, and delta
}
\description{
MVP.BRENT.Vg.Ve
variance component estimation using the BRENT method
}
\examples{
\donttest{
phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
phenotype <- read.table(phePath, header=TRUE)
print(dim(phenotype))
genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
genotype <- attach.big.matrix(genoPath)
print(dim(genotype))

eigenK <- eigen(MVP.K.VanRaden(genotype, cpu=1))
vc <- MVP.BRENT.Vg.Ve(y=phenotype[,2], X=matrix(1, nrow(phenotype)), eigenK=eigenK)
print(vc)
}

}
