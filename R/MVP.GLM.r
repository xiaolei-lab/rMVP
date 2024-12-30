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
#' Last update: May 25, 2017
#'
#' @author Lilin Yin and Xiaolei Liu
#' 
#' @param phe phenotype, n * 2 matrix
#' @param geno genotype, either m by n or n by m is supportable, m is marker size, n is population size
#' @param CV Covariance, design matrix(n * x) for the fixed effects
#' @param ind_idx the index of effective genotyped individuals
#' @param mrk_idx the index of effective markers used in analysis
#' @param mrk_bycol whether the markers are stored by columns in genotype (i.e. M is a n by m matrix)
#' @param maxLine the number of markers handled at a time, smaller value would reduce the memory cost
#' @param cpu number of cpus used for parallel computation
#' @param verbose whether to print detail.
#'
#' @return m * 2 matrix, the first column is the SNP effect, the second column is the P values
#' @export
#'
#' @examples
#' \donttest{
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' idx <- !is.na(phenotype[, 2])
#' phenotype <- phenotype[idx, ]
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' genotype <- deepcopy(genotype, rows=idx)
#' print(dim(genotype))
#' 
#' glm <- MVP.GLM(phe=phenotype, geno=genotype, cpu=1)
#' str(glm)
#' }
MVP.GLM <-
function(
    phe, 
    geno, 
    CV=NULL, 
    ind_idx=NULL,
    mrk_idx=NULL,
    mrk_bycol=TRUE,
    maxLine=5000,
    cpu=1,
    verbose=TRUE
){

    if(is.null(ind_idx)){
        if(nrow(phe) != ncol(geno) && nrow(phe) != nrow(geno)) stop("number of individuals does not match in phenotype and genotype.")
        n <- ifelse(nrow(phe) == ncol(geno), ncol(geno), nrow(geno))
        if(mrk_bycol){
            if(nrow(phe) != nrow(geno))    stop("set the argument 'mrk_bycol=FALSE'?")
        }else{
            if(nrow(phe) != ncol(geno))    stop("set the argument 'mrk_bycol=TRUE'?")
        }
    }else{
        n <- length(ind_idx)
        if(nrow(phe) != n) stop("number of individuals does not match in phenotype and genotype.")
    }

    if(!is.big.matrix(geno))    stop("genotype should be in 'big.matrix' format.")
    ys <- as.numeric(as.matrix(phe[,2]))
    if(sum(is.na(ys)) != 0) stop("NAs are not allowed in phenotype.")
    
    if(is.null(CV)){
        X0 <- matrix(1, n)
    }else{
        if(nrow(CV) != n)   stop("number of individuals does not match in phenotype and fixed effects.")
        if(sum(is.na(CV)) != 0) stop("NAs are not allowed in fixed effects.")
        CV.index <- apply(CV, 2, function(x) length(table(x)) > 1)
        CV <- CV[, CV.index, drop=FALSE]
        X0 <- cbind(matrix(1, n), CV)
    }
    X0 <- as.matrix(X0)
    iX0X0 <- MASS::ginv(crossprod(X0))

    logging.log("scanning...\n", verbose = verbose)

    mkl_env({
        results <- glm_c(y = ys, X = X0, iXX = iX0X0, geno@address, geno_ind = ind_idx, marker_ind = mrk_idx, marker_bycol = mrk_bycol, step = maxLine, verbose = verbose, threads = cpu)
    })

    return(results[, c(1, 2, ncol(results))])
}#end of MVP.GLM function
