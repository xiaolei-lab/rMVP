# Data pre-processing module
# 
# Copyright (C) 2016-2018 by Xiaolei Lab
# 
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


#' To perform GWAS with GLM and MLM model and get the P value of SNPs
#'
#' Build date: Aug 30, 2016
#' Last update: Aug 30, 2016
#' 
#' @author Lilin Yin and Xiaolei Liu
#' 
#' @param phe phenotype, n * 2 matrix
#' @param geno genotype, m * n, m is marker size, n is population size
#' @param K Kinship, Covariance matrix(n * n) for random effects; must be positive semi-definite
#' @param CV covariates
#' @param REML a list that contains ve and vg
#' @param priority speed or memory
#' @param cpu number of cpus used for parallel computation
#' @param bar whether to show the progress bar
#' @param vc.method the methods for estimating variance component("emma" or "gemma")
#' @param maxLine when the priority is 'memory', users can change this parameter to limit the memory
#' @param file.output whether to output files or not
#' @param memo a marker on output file name
#'
#' @return
#' results: a m * 2 matrix, the first column is the SNP effect, the second column is the P values
#' @export
#'
#' @examples
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' idx <- !is.na(phenotype[, 2])
#' phenotype <- phenotype[idx, ]
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' genotype <- genotype[, idx]
#' print(dim(genotype))
#' K <- MVP.K.VanRaden(genotype)
#' mlm <- MVP.MLM(phe=phenotype, geno=genotype, K=K)
#' str(mlm)
MVP.MLM <-
function(phe, geno, K=NULL, CV=NULL, REML=NULL, priority="speed", cpu=1, bar=TRUE,vc.method="EMMA",maxLine=1000, file.output=TRUE, memo="MVP"){
    if (Sys.info()[['sysname']] == 'Windows')
        cpu <- 1
    math.cpu <- try(getMKLthreads(), silent=TRUE)
    
    n <- ncol(geno)
    m <- nrow(geno)
    if(priority=="speed") geno <- as.matrix(geno)
    
    ys <- as.numeric(as.matrix(phe[,2]))
    if (is.null(K)) {
        print("Calculating Kinship...")
        K <- MVP.K.VanRaden(M=geno, priority=priority, maxLine=maxLine);gc()
        if (file.output) {
            filebck <- paste0("MVP.", colnames(phe)[2], memo, ".kin.bin")
            filedes <- paste0("MVP.", colnames(phe)[2], memo, ".kin.desc")
            if (file.exists(filebck)) file.remove(filebck)
            if (file.exists(filedes)) file.remove(filedes)
            Kin.backed <- big.matrix(
                nrow = nrow(K),
                ncol = ncol(K),
                type = "double",
                backingfile = filebck,
                descriptorfile = filedes
            )
            Kin.backed[, ] <- K[, ]
            flush(Kin.backed)
            rm(list = c("Kin.backed"))
            gc()
        }
    } else {
        # convert K to base:matrix
        K <- K[, ]
    }

    if (is.null(CV)) {
        X0 <- matrix(1, n)
    } else {
        X0 <- cbind(matrix(1, n), CV)
    }

    # number of fixed effects
    nf <- ncol(X0) + 1
    if (is.null(REML)) {
        print("Variance components...")   
        if (vc.method == "EMMA")
            REML <- MVP.EMMA.Vg.Ve(y=ys, X=X0, K=K)
        if (vc.method == "GEMMA")
            REML <- MVP.GEMMA.Vg.Ve(y=ys, X=X0, K=K)
        print("Variance components is Done!")
    }

    q0 <- ncol(X0)
    iXX <- matrix(NA,q0+1,q0+1)
    Xt <- matrix(NA,n, q0+1)
 
    ves <- REML$ve
    vgs <- REML$vg
    lambda <- ves/vgs
    print("Eigen Decomposition...")
    eig <- eigen(K, symmetric=TRUE); rm(K); gc()
    print("Eigen-Decomposition is Done!")
    
    U <- eig$vectors * matrix(sqrt(1/(eig$values + lambda)), n, length(eig$values), byrow=TRUE); rm(eig); gc()
    
    y <- matrix(ys)
    yt <- crossprod(U, y)
    X0t <- crossprod(U, X0)
    X0X0 <- crossprod(X0t)
    X0Y <- crossprod(X0t,yt)
    iX0X0 <- solve(X0X0)
    
    Xt[1:n,1:q0] <- X0t

    #parallel function for MLM model
    eff.mlm.parallel <- function(i){
        if(bar) print.f(i)
        # if(i%%1000==0){
            # print(paste("****************", i, "****************",sep=""))
        # }
        #if(cpu>1 & r.open) setMKLthreads(math.cpu)

        SNP <- geno[i, ]
        xst <- crossprod(U, SNP)
        Xt[1:n,q0+1] <- xst
        X0Xst <- crossprod(X0t,xst)
        XstX0 <- t(X0Xst)
        xstxst <- crossprod(xst)
        xsY <- crossprod(xst,yt)
        XY <- c(X0Y,xsY)
        #B22 <- xstxst - XstX0%*%iX0X0%*%X0Xst
        B22 <- xstxst - crossprod(X0Xst, iX0X0) %*% X0Xst
        invB22 <- 1/B22
        B21 <- tcrossprod(XstX0, iX0X0)
        NeginvB22B21 <- crossprod(-invB22,B21)
        B11 <- iX0X0 + as.numeric(invB22)*crossprod(B21,B21)
            
        iXX[1:q0,1:q0]=B11
        iXX[(q0+1),(q0+1)]=1/B22
        iXX[(q0+1),1:q0]=NeginvB22B21
        iXX[1:q0,(q0+1)]=NeginvB22B21
        beta <- crossprod(iXX,XY)
        stats <- beta[(q0+1)]/sqrt((iXX[(q0+1), (q0+1)]) * vgs)
        p <- 2 * pt(abs(stats), n-(q0+1), lower.tail=FALSE)
        effect<- beta[(q0+1)]
        return(list(effect = effect, p = p))
    }
    
    #Paralleled MLM
    if (cpu == 1) {
        print.f <- function(i) {
            print_bar(i = i, n = m, type = "type1", fixed.points = TRUE)
        }
        mkl_env({
            results <- lapply(seq_len(m), eff.mlm.parallel)
        }, threads = min(2^(n %/% 1000), math.cpu))
    } else {
        if (Sys.info()[['sysname']] == 'Windows'){
            print.f <- function(i) {
                print_bar(i = i, n = m, type = "type1", fixed.points = TRUE)
            }
            cl <- makeCluster(getOption("cl.cores", cpu))
            clusterExport(cl, varlist=c("geno", "yt", "X0", "U", "vgs", "ves", "math.cpu"), envir=environment())
            Exp.packages <- clusterEvalQ(cl, c(library(bigmemory)))
            results <- parLapply(cl, 1:m, eff.mlm.parallel)
            stopCluster(cl)
        }else{
            tmpf.name <- tempfile()
            tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
            writeBin(0, tmpf)
            print.f <- function(i){print_bar(n=m, type="type3", tmp.file=tmpf, fixed.points=TRUE)}
            R.ver <- Sys.info()[['sysname']]
            if(R.ver == 'Linux') {
                math.cpu <- try(getMKLthreads(), silent=TRUE)
                try(setMKLthreads(1), silent=TRUE)
            }
            results <- mclapply(1:m, eff.mlm.parallel, mc.cores=cpu)
            if(R.ver == 'Linux') {
                try(setMKLthreads(math.cpu), silent=TRUE)
                #try(setMKLthreads(1), silent=TRUE)
            }
            close(tmpf); unlink(tmpf.name); cat('\n');
        }
    }
    if(is.list(results)) results <- matrix(unlist(results), m, byrow=TRUE)
    #print("****************MLM ACCOMPLISHED****************")
    return(results)
}#end of MVP.MLM function
