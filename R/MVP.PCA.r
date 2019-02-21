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


#' Principal Component Analysis
#' 
#' Build date: Dec 14, 2016
#' Last update: Oct 29, 2018
#' @author Xiaolei Liu, Lilin Yin and Haohao Zhang
#' 
#' @param M Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param perc percentage of total SNPs selected for calculate principal components
#' @param pcs.keep maximum number of PCs for output
#' @param memo a marker on output file name
#'
#' @return
#' Output: PCs - a n * npc matrix of top number of PCs, n is population size and npc is @param pcs.keep 
#' 
#' @export
#'
#' @examples
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' pca <- MVP.PCA(M=genotype)
#' str(pca)
MVP.PCA <-
function(M, perc=1, pcs.keep=5, memo=NULL){
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    linux <- R.ver == 'Linux'
    mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent = TRUE),"try-error")

    if (r.open && mac) {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    #Data Check
    if (is.null(M)) {
        stop("There is no genotype data!")
    }

    if (pcs.keep < 0) {
        pcs.keep = 0
    }

    if (pcs.keep == 0) {
        return(NULL)
    }

    if (perc < 0 | perc > 1) {
        print("perc value should less than 1 and bigger than 0, MVP reset it to 1")
        perc <- 1
    }
    print("Principal Component Analysis Start...")
    if (perc < 1) {
        # set.seed(123456)
        m <- nrow(M)
        big.geno <- t(M[sample(1:m, round(m * perc)), ])
    }else{
        big.geno <- t(as.matrix(M))
    }
    #Calculate PCs by using bigPCA package

    PCs <- prcomp(big.geno)$x[, 1:pcs.keep]

    return(list(PCs=PCs))
} # end of MVP.PCA function
