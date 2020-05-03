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
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param CV Covariance, design matrix(n * x) for the fixed effects
#' @param cpu number of cpus used for parallel computation
#' @param verbose whether to print detail.
#'
#' @return m * 2 matrix, the first column is the SNP effect, the second column is the P values
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
#' genotype <- deepcopy(genotype, cols=idx)
#' print(dim(genotype))
#' 
#' glm <- MVP.GLM(phe=phenotype, geno=genotype)
#' str(glm)
MVP.GLM <-
function(
    phe, 
    geno, 
    CV=NULL, 
    cpu=1, 
    verbose=TRUE
){

    n <- ncol(geno)
    m <- nrow(geno)
    ys <- as.numeric(as.matrix(phe[,2]))
    
    if(!is.big.matrix(geno))    stop("genotype should be in 'big.matrix' format.")
    if(sum(is.na(ys)) != 0) stop("NAs are not allowed in phenotype.")
    if(nrow(phe) != ncol(geno)) stop("number of individuals not match in phenotype and genotype.")
    
    if(is.null(CV)){
        X0 <- matrix(1, n)
    }else{
        CV.index <- apply(CV, 2, function(x) length(table(x)) > 1)
        CV <- CV[, CV.index, drop=FALSE]
        X0 <- cbind(matrix(1, n), CV)
    }
    X0 <- as.matrix(X0)
    iX0X0 <- MASS::ginv(crossprod(X0))

    logging.log("scanning...\n", verbose = verbose)

    mkl_env({
        results <- glm_c(y = ys, X = X0, iXX = iX0X0, geno@address, verbose = verbose, threads = cpu)
    })

    return(results[, c(1, 2, ncol(results))])
}#end of MVP.GLM function
