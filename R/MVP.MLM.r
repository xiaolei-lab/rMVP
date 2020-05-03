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
#' idx <- !is.na(phenotype[, 2])
#' phenotype <- phenotype[idx, ]
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' genotype <- deepcopy(genotype, cols=idx)
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
    vc.method=c("BRENT", "EMMA", "HE"),
    verbose=TRUE
){
    R.ver <- Sys.info()[['sysname']]
    r.open <- eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))
    
    if (R.ver == 'Windows') cpu <- 1
    if (r.open && cpu > 1 && R.ver == 'Darwin') {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }

    vc.method <- match.arg(vc.method)
    n <- ncol(geno)
    m <- nrow(geno)
    ys <- as.numeric(as.matrix(phe[,2]))

    if(!is.big.matrix(geno))    stop("genotype should be in 'big.matrix' format.")
    if(sum(is.na(ys)) != 0) stop("NAs are not allowed in phenotype.")
    if(nrow(phe) != ncol(geno)) stop("number of individuals not match in phenotype and genotype.")
    
    if(is.null(K)){
        if(vc.method == "EMMA" | vc.method == "he")  stop("Kinship must be provided!")
        if(vc.method == "BRENT"){
            if(is.null(eigenK)) stop("eigenK must be provided!")
        }
    }else{
        # convert K to base:matrix
        K <- K[, ]
        if(is.null(eigenK)){
            logging.log("Eigen Decomposition on Genomic Relationship Matrix", "\n", verbose = verbose)
            eigenK <- eigen(K, symmetric=TRUE)
        }
    }

    if (is.null(CV)) {
        X0 <- matrix(1, n)
    } else {
        CV.index <- apply(CV, 2, function(x) length(table(x)) > 1)
        CV <- CV[, CV.index, drop=FALSE]
        X0 <- cbind(matrix(1, n), CV)
    }
    X0 <- as.matrix(X0)

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

    ves <- REML$ve
    vgs <- REML$vg
    lambda <- ves/vgs

    U <- eigenK$vectors * matrix(sqrt(1/(eigenK$values + lambda)), n, length(eigenK$values), byrow=TRUE); rm(eigenK); gc()
    
    logging.log("scanning...\n", verbose = verbose)
    mkl_env({
        results <- mlm_c(y = ys, X = X0, U = U, vgs = vgs, geno@address, verbose = verbose, threads = cpu)
    })

    return(results)
}#end of MVP.MLM function
