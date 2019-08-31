
#' MVP.BRENT.Vg.Ve
#' variance component estimation using the BRENT method
#' 
#' @param y phenotype
#' @param X covariate matrix, the first column is 1s
#' @param eigenK eigen of Kinship matrix
#' @param verbose whether to print detail.
#' 
#' @return vg, ve, and delta
#' @export
#'
#' @examples
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' 
#' eigenK <- eigen(MVP.K.VanRaden(genotype))
#' vc <- MVP.BRENT.Vg.Ve(y=phenotype[,2], X=matrix(1, nrow(phenotype)), eigenK=eigenK)
#' print(vc)
#' 
MVP.BRENT.Vg.Ve <- function(y, X, eigenK, verbose = FALSE) {
    p = 0
    Sigma <- eigenK$values
    w <- which(Sigma < 1e-6)
    Sigma[w] <- 1e-6
    U <- eigenK$vectors
    min_h2 = 0
    max_h2 = 1
    tol = .Machine$double.eps^0.25
    reml <- fit_diago_brent(y, X, p, Sigma, U, min_h2, max_h2, tol, verbose = verbose)
    vg <- reml[[2]]
    ve <- reml[[1]]
    delta <- ve / vg
    return(list(vg = vg, ve = ve, delta = delta))
}
