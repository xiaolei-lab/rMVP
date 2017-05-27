MVP.EMMA.Vg.Ve <-
function(y, X, K, ngrids=100, llim=-10, ulim=10, esp=1e-10){
##########################################################################################################
# Object: Estimate variance components using EMMA
# Input:
# y: phenotype
# X: covariate matrix, the first column is 1s
# K: Kinship matrix
# ngrids, llim, ulim, and esp: parameters for estimating vg and ve
# Authors: EMMA (Kang et. al. Genetics, 2008)
# Modified only for speed up by Xiaolei Liu and Lilin Yin
# Build date: August 30, 2016
# Last update: January 27, 2017
##########################################################################################################

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
            iXX <- geninv(XX)
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
    print("Variance Components Estimation is Done!")
    return (list(REML=maxLL, delta=maxdelta, ve=maxve, vg=maxva))
}# end of MVP.EMMA.Vg.Ve function
