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
#' @param eigenK list of eigen Kinship
#' @param CV covariates
#' @param REML a list that contains ve and vg
#' @param cpu number of cpus used for parallel computation
#' @param bar whether to show the progress bar
#' @param vc.method the methods for estimating variance component("emma" or "he" or "brent")
#' @param verbose whether to print detail.
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
#' 
#' mlm <- MVP.MLM(phe=phenotype, geno=genotype, K=K)
#' str(mlm)
#' 

MVP.MLM <-
function(
    phe, 
    geno, 
    K=NULL,
    eigenK=NULL,
    CV=NULL, 
    REML=NULL,
    cpu=1,
    bar=TRUE,
    vc.method=c("BRENT", "EMMA", "HE"),
    verbose=TRUE
){
    R.ver <- Sys.info()[['sysname']]
    r.open <- eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))
    
    if (R.ver == 'Windows') cpu <- 1
    if (r.open && cpu > 1 && R.ver == 'Darwin') {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }

    math.cpu <- eval(parse(text = "try(getMKLthreads(), silent=TRUE)"))
    vc.method <- match.arg(vc.method)
    n <- ncol(geno)
    m <- nrow(geno)
    
    ys <- as.numeric(as.matrix(phe[,2]))
    if(is.null(K)){
        if(vc.method == "EMMA" | vc.method == "he")  stop("Kinship must be provided!")
        if(vc.method == "BRENT"){
            if(is.null(eigenK)) stop("eigenK must be provided!")
        }
    }else{
        # convert K to base:matrix
        K <- K[, ]
        if(is.null(eigenK)){
            logging.log("Eigen Decomposition of Kinship...", "\n", verbose = verbose)
            eigenK <- eigen(K, symmetric=TRUE)
            logging.log("Eigen-Decomposition is Done!", "\n", verbose = verbose)
        }
    }

    if (is.null(CV)) {
        X0 <- matrix(1, n)
    } else {
	CV.index <- apply(CV, 2, function(x) length(table(x)) > 1)
	CV <- CV[, CV.index, drop=FALSE]
        X0 <- cbind(matrix(1, n), CV)
    }

    # number of fixed effects
    nf <- ncol(X0) + 1
    if(is.null(REML)) {
		logging.log(paste("Variance components using: ", vc.method, sep=""), "\n", verbose = verbose)   
        if (vc.method == "EMMA") REML <- MVP.EMMA.Vg.Ve(y=ys, X=X0, K=K)
        if (vc.method == "HE") REML <- MVP.HE.Vg.Ve(y=ys, X=X0, K=K)
        if (vc.method == "BRENT") REML <- MVP.BRENT.Vg.Ve(y=ys, X=X0, eigenK=eigenK)
		logging.log(paste("Estimated Vg and Ve: ", sprintf("%.6f", REML$vg), " ", sprintf("%.6f", REML$ve), sep=""), "\n", verbose = verbose)
    }else{
        logging.log(paste("Provided Vg and Ve: ", sprintf("%.6f", REML$vg), " ", sprintf("%.6f", REML$ve), sep=""), "\n", verbose = verbose)
    }
    if(!is.null(K)){rm(K); gc()}

    q0 <- ncol(X0)
    iXX <- matrix(NA,q0+1,q0+1)
    Xt <- matrix(NA,n, q0+1)
 
    ves <- REML$ve
    vgs <- REML$vg
    lambda <- ves/vgs

    U <- eigenK$vectors * matrix(sqrt(1/(eigenK$values + lambda)), n, length(eigenK$values), byrow=TRUE); rm(eigenK); gc()
    
    y <- matrix(ys)
    yt <- crossprod(U, y)
    X0t <- crossprod(U, X0)
    X0X0 <- crossprod(X0t)
    X0Y <- crossprod(X0t,yt)
    iX0X0 <- ginv(X0X0)
    #iX0X0 <- solve(X0X0)    
    Xt[1:n,1:q0] <- X0t

    #parallel function for MLM model
    eff.mlm.parallel <- function(i){
        if(bar) print.f(i)
        SNP <- geno[i, ]
        xst <- crossprod(U, SNP)
        Xt[1:n,q0+1] <- xst
        X0Xst <- crossprod(X0t,xst)
        XstX0 <- t(X0Xst)
        xstxst <- crossprod(xst)
        xsY <- crossprod(xst,yt)
        XY <- c(X0Y,xsY)
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
        se <- sqrt((iXX[(q0+1), (q0+1)]) * vgs)
        stats <- beta[(q0+1)]/se
        p <- 2 * pt(abs(stats), n-(q0+1), lower.tail=FALSE)
        effect<- beta[(q0+1)]
        return(list(effect = effect, se = se, p = p))
    }
    logging.log("scanning...\n", verbose = verbose)
    #Paralleled MLM
    if(cpu == 1){
        math.cpu <- eval(parse(text = "try(getMKLthreads(), silent=TRUE)"))
        mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
        eval(parse(text = "try(setMKLthreads(mkl.cpu), silent=TRUE)"))
        print.f <- function(i){print_bar(i=i, n=m, type="type1", fixed.points=TRUE, verbose = verbose)}
        results <- lapply(1:m, eff.mlm.parallel)
        eval(parse(text = "try(setMKLthreads(math.cpu), silent=TRUE)"))
    }else{
        if(R.ver == 'Windows'){
            print.f <- function(i){print_bar(i=i, n=m, type="type1", fixed.points=TRUE, verbose = verbose)}
            cl <- makeCluster(getOption("cl.cores", cpu))
            clusterExport(cl, varlist=c("geno", "yt", "X0", "U", "vgs", "ves", "math.cpu"), envir=environment())
            Exp.packages <- clusterEvalQ(cl, c(library(bigmemory)))
            results <- parLapply(cl, 1:m, eff.mlm.parallel)
            stopCluster(cl)
        }else{
            tmpf.name <- tempfile()
            tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
            writeBin(0, tmpf)
            print.f <- function(i){print_bar(n=m, type="type3", tmp.file=tmpf, fixed.points=TRUE, verbose = verbose)}
            mkl_env({
                results <- mclapply(1:m, eff.mlm.parallel, mc.cores=cpu)
            })
            close(tmpf); unlink(tmpf.name); logging.log('\n', verbose = verbose)
        }
    }
    if(is.list(results)) results <- matrix(unlist(results), m, byrow=TRUE)
    return(results)
}#end of MVP.MLM function

