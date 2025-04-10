% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MVP.Data.r
\name{MVP.Data.Kin}
\alias{MVP.Data.Kin}
\title{Kinship}
\usage{
MVP.Data.Kin(
  fileKin = TRUE,
  mvp_prefix = "mvp",
  out = NULL,
  maxLine = 10000,
  mrk_bycol = TRUE,
  sep = "\\t",
  cpu = 1,
  verbose = TRUE
)
}
\arguments{
\item{fileKin}{Kinship that represents relationship among individuals, n * n matrix, n is sample size}

\item{mvp_prefix}{Prefix for mvp format files}

\item{out}{prefix of output file name}

\item{maxLine}{the number of markers handled at a time, smaller value would reduce the memory cost}

\item{mrk_bycol}{whether the markers are stored by columns in genotype (i.e. genotype is a n by m matrix)}

\item{sep}{seperator for Kinship file.}

\item{cpu}{the number of cpu}

\item{verbose}{whether to print detail.}
}
\value{
Output file:
<out>.kin.bin
<out>.kin.desc
}
\description{
Kinship
}
\examples{
\donttest{
geno <- file.path(system.file("extdata", "06_mvp-impute", package = "rMVP"), "mvp.imp")

MVP.Data.Kin(TRUE, mvp_prefix=geno, out=tempfile("outfile"), cpu=1)
}

}
