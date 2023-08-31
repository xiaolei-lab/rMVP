# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Estimate variance components using EMMA
#'
#' Build date: August 30, 2016
#' Last update: January 27, 2017
#' 
#' @author EMMA (Kang et. al. Genetics, 2008), Modified only for speed up by Xiaolei Liu and Lilin Yin
#' 
#' @param y phenotype, n * 2
#' @param X covariate matrix, the first column is 1s
#' @param K Kinship matrix
#' @param ngrids parameters for estimating vg and ve
#' @param llim parameters for estimating vg and ve
#' @param ulim parameters for estimating vg and ve
#' @param esp parameters for estimating vg and ve
#'
#' @return
#' Output: REML - maximum log likelihood
#' Output: delta - exp(root)
#' Output: ve - residual variance
#' Output: vg - genetic variance
#' 
#' @export
#'
#' @examples
#' \donttest{
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' 
#' K <- MVP.K.VanRaden(genotype, cpu=1)
#' vc <- MVP.EMMA.Vg.Ve(y=phenotype[,2], X=matrix(1, nrow(phenotype)), K=K)
#' print(vc)
#' }
#' 
MVP.EMMA.Vg.Ve <-
function(y, X, K, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
    # NA in phenotype
    idx <- !is.na(y)
    y <- y[idx]
    X <- as.matrix(X[idx, ])
    K <- K[idx, idx]
    
    if(!is.numeric(y))	y <- as.numeric(as.character(y))
    emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas)
    {
        nq <- length(etas)
        delta <-  exp(logdelta)
        return( 0.5 * (nq * (log(nq/(2 * pi))-1-log(sum(etas * etas/(lambda + delta))))-sum(log(lambda + delta))) )
    }
    emma.eigen.R.wo.Z=function(K, X) {
        n <- nrow(X)
        q <- ncol(X)
        
        XX <- crossprod(X)
        iXX <- try(solve(XX), silent = TRUE)
        if(inherits(iXX, "try-error")){
            #library(MASS)
            iXX <- ginv(XX)
        }

        SS1 <- X %*% iXX
        SS2 <- tcrossprod(SS1, X)
        S <- diag(n)-SS2
        
        eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric=TRUE)#S4 error here

        stopifnot(!is.complex(eig$values))
        return(list(values=eig$values[1:(n-q)]-1, vectors=eig$vectors[, 1:(n-q)]))
    }
    eig.R=emma.eigen.R.wo.Z(K=K, X=X)
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)

    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)

    etas <- crossprod(eig.R$vectors, y)
    
    logdelta <- (0:ngrids)/ngrids * (ulim-llim) + llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values, n-q, m) + matrix(delta, n-q, m, byrow=TRUE)
    Etasq <- matrix(etas * etas, n-q, m)
    LL <- 0.5 * ((n-q) * (log((n-q)/(2 * pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5 * delta * ((n-q) * colSums(Etasq/(Lambdas * Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim, eig.R$values, etas))
    }
    if( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim, eig.R$values, etas))
    }
    for(i in 1:(m-1) ){
        if( ( dLL[i] * dLL[i + 1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i + 1] < 0 ) ) {
            emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
                nq <- length(etas)
                delta <- exp(logdelta)
                etasq <- etas * etas
                ldelta <- lambda + delta
                return( 0.5 * (nq * sum(etasq/(ldelta * ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
            }
            r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i + 1], lambda=eig.R$values, etas=etas)
            optlogdelta <- append(optlogdelta, r$root)
            emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
                nq <- length(etas)
                delta <-  exp(logdelta)
                return( 0.5 * (nq * (log(nq/(2 * pi))-1-log(sum(etas * etas/(lambda + delta))))-sum(log(lambda + delta))) )
            }
            optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root, eig.R$values, etas))
        }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    #handler of grids with NaN log
    replaceNaN<-function(LL) {
        index=(LL == "NaN")
        if(length(index)>0) theMin=min(LL[!index])
        if(length(index)<1) theMin="NaN"
        LL[index]=theMin
        return(LL)    
    }
    optLL=replaceNaN(optLL)   
    maxLL <- max(optLL)
    maxva <- sum(etas * etas/(eig.R$values + maxdelta))/(n-q)    
    maxve <- maxva * maxdelta
    return (list(REML=maxLL, delta=maxdelta, ve=maxve, vg=maxva))
}# end of MVP.EMMA.Vg.Ve function
