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
#' @param K kinship matrix
#' @param priority speed or memory
#' @param pcs.keep maximum number of PCs for output
#' @param cpu the number of cpu
#' @param verbose whether to print detail.
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
#' 
#' pca <- MVP.PCA(M=genotype)
#' str(pca)
#' 
MVP.PCA <-
function(M=NULL, K=NULL, priority=c("speed", "memory"), pcs.keep=5, cpu=1, verbose=TRUE){
    # R.ver <- Sys.info()[['sysname']]
    # wind <- R.ver == 'Windows'
    # linux <- R.ver == 'Linux'
    # mac <- (!linux) & (!wind)
    # r.open <- eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))

    # if (r.open && mac) {
    #     Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    # }
    #Data Check
    if (is.null(M) & is.null(K)) {
        stop("There is no genotype data or relationship matrix!")
    }

    priority <- match.arg(priority)

    if (pcs.keep < 0) {
        pcs.keep = 0
    }

    if (pcs.keep == 0) {
        return(NULL)
    }

    if(is.null(K)){
        K <- MVP.K.VanRaden(M=M, priority=priority, cpu=cpu)
    }

    logging.log("Eigen Decomposition on GRM", "\n", verbose = verbose)
    # if(r.open)  eval(parse(text = "try(setMKLthreads(cpu), silent=TRUE)"))
    PCs <- eigen(K, symmetric=TRUE)$vectors[, 1:pcs.keep]
    logging.log("Deriving PCs successfully", "\n", verbose = verbose)

    return(PCs=PCs)
}#end of MVP.PCA function

