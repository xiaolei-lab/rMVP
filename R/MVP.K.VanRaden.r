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


#' Calculate Kinship matrix by VanRaden method
#'
#' @param M genotype, either m by n or n by m is supportable, m is marker size, n is population size
#' @param maxLine the number of markers handled at a time, smaller value would reduce the memory cost
#' @param ind_idx the index of effective genotyped individuals used in analysis
#' @param mrk_idx the index of effective markers used in analysis
#' @param mrk_freq the prior calculated major allele frequency (not MAF) for all markers used in analysis
#' @param mrk_bycol whether the markers are stored by columns in genotype (i.e. M is a n by m matrix)
#' @param cpu the number of cpu
#' @param verbose whether to print detail.
#' @param checkNA whether to check NA in genotype.
#'
#' @return K, n * n matrix
#' @export
#'
#' @examples
#' \donttest{
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#'
#' K <- MVP.K.VanRaden(genotype, cpu=1)
#' }
#'
MVP.K.VanRaden <-
function(
    M,
    maxLine=5000,
    ind_idx=NULL,
    mrk_idx=NULL,
    mrk_freq = NULL,
    mrk_bycol = TRUE,
    cpu=1,
    verbose=TRUE,
    checkNA=TRUE
){
    
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    linux <- R.ver == 'Linux'
    mac <- (!linux) & (!wind)
    HPCMathLib <- mac | grepl("mkl", sessionInfo()$LAPACK) | grepl("openblas", sessionInfo()$LAPACK) | eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))

    if(!is.big.matrix(M)) stop("Format of Genotype Data must be big.matrix")
    if(checkNA){
        if(hasNA(M@address, mrkbycol = mrk_bycol, geno_ind = ind_idx, marker_ind = mrk_idx, threads = cpu))   stop("NA is not allowed in genotype, use 'MVP.Data.impute' to impute.")
    }

    if(!is.null(mrk_freq)){
        if(!is.null(mrk_idx)){
            if(length(mrk_freq) != length(mrk_idx)) stop("mismatched number of markers between provided frequency and index.")
        }else{
            if(mrk_bycol){
                if(length(mrk_freq) != ncol(M)) stop("mismatched number of markers between genotype and provided frequency.")
            }else{
                if(length(mrk_freq) != nrow(M)) stop("mismatched number of markers between genotype and provided frequency.")
            }
        }
    }
    
    n <- ifelse(is.null(ind_idx), ifelse(mrk_bycol, nrow(M), ncol(M)), length(ind_idx))
    m <- ifelse(is.null(mrk_idx), ifelse(mrk_bycol, ncol(M), nrow(M)), length(mrk_idx))

    # logging.log("Relationship matrix mode in", priority[1], "\n", verbose = verbose)
    # if(is.null(dim(M))) M <- t(as.matrix(M))
    
    logging.log(paste0("Computing GRM for ", n, " individuals using ", m, " markers with a step of ", maxLine), "\n", verbose = verbose)
    
    K <- try(kin_cal(M@address, geno_ind = ind_idx, marker_ind = mrk_idx, marker_freq = mrk_freq, marker_bycol = mrk_bycol, threads = cpu, step = maxLine, verbose = verbose, mkl = HPCMathLib), silent=TRUE)
    # K <- try(kin_cal_s(M@address, threads = cpu, verbose = verbose, mkl = HPCMathLib), silent=TRUE)

    if(inherits(K,"try-error")){
        logging.log("Out of memory, please set a smaller value to parameter 'maxLine' and try again.", "\n", verbose = verbose)
        stop(K[[1]])
    }
    
    # switch(
    #     match.arg(priority),
    #     "speed" = {
    #         # if (!is.matrix(M)) M <- as.matrix(M)
    #         # n <- ncol(M)
    #         # m <- nrow(M)
    #         # Pi <- 0.5 * rowMeans(M)
    #         # logging.log("Scale the genotype matrix", "\n", verbose = verbose)
    #         # M <- M - 2 * Pi
    #         # SUM <- sum(Pi * (1 - Pi))


    #         K <- try(kin_cal_s(M@address, threads = cpu, verbose = verbose, mkl = HPCMathLib), silent=TRUE)

    #         if(inherits(K,"try-error")){
    #             logging.log("Out of memory, please set parameter (..., priority='memory') and try again.", "\n", verbose = verbose)
    #             stop(K[[1]])
    #         }
    #     },

    #     "memory" = {
    #         K <- kin_cal_m(M@address, threads=cpu, verbose = verbose)
    #         # n <- ncol(M)
    #         # m <- nrow(M)
    #         # bac <- paste0("Z", memo, ".temp.bin")
    #         # des <- paste0("Z", memo, ".temp.desc")
    #         # if (file.exists(bac)) file.remove(bac)
    #         # if (file.exists(des)) file.remove(des)
    #         # #options(bigmemory.typecast.warning=FALSE)
    #         # Z <- big.matrix(
    #         #     nrow = m,
    #         #     ncol = n,
    #         #     type = "double",
    #         #     backingfile = bac,
    #         #     descriptorfile = des,
    #         #     init = 0.1
    #         # )
    #         # Pi <- NULL
    #         # estimate.memory <- function(dat, integer=FALSE, raw=FALSE){
    #         #     cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
    #         #     dimz <- dat
    #         #     if(length(dimz) == 1) { dimz[2] <- 1 }
    #         #     if(length(dimz)>1 & length(dimz)<11 & is.numeric(dimz)) {
    #         #         total.size <- as.double(1)
    #         #         for(cc in 1:length(dimz)) { total.size <- as.double(total.size * as.double(dimz[cc])) }
    #         #         memory.estimate <- as.double(as.double(total.size)/cells.per.gb)
    #         #         memory.estimate <- memory.estimate
    #         #         if(integer) { memory.estimate <- memory.estimate/2 } else { if(raw) { memory.estimate <- memory.estimate/8 } }
    #         #         return(memory.estimate)
    #         #     } else {
    #         #         # guessing this is a vector
    #         #         if(!is.list(dimz) & is.vector(dimz)) {
    #         #             LL <- length(dimz)
    #         #             return(estimate.memory(LL, integer=integer, raw=raw))
    #         #         } else {
    #         #             warning("tried to estimate memory for object which is neither a vector, pair of dimension sizes or a dataframe/matrix")
    #         #         }
    #         #     }
    #         # }
    #         # if((Sys.info()[['sysname']]) == 'Windows'){
    #         #     max.gb <- memory.limit()/1000
    #         # }else{
    #         #     max.gb <- Inf
    #         # }
    #         # maxLines.gb <- estimate.memory(c(maxLine, n))
    #         # if(maxLines.gb > max.gb) stop("Memory limited! Please reset the 'maxLine'")
    #         # loop.index <- seq(0, m, maxLine)[-1]
    #         # if(max(loop.index) < m) loop.index <- c(loop.index, m)
    #         # loop.len <- length(loop.index)
    #         # print("Z assignment...")
    #         # for(cc in 1:loop.len){
    #         #     if(loop.len == 1){
    #         #         c1 <- 1
    #         #     }else{
    #         #         c1 <- ifelse(cc == loop.len, (loop.index[cc-1]) + 1, loop.index[cc]-maxLine + 1)
    #         #     }
    #         #     c2 <- loop.index[cc]
    #         #     means <-rowMeans(M[c1:c2, 1:n])
    #         #     if(!is.null(weight)){
    #         #         Z[c1:c2, 1:n] <- (M[c1:c2, 1:n]-means) * sqrt(weight[c1:c2])
    #         #     }else{
    #         #         Z[c1:c2, 1:n] <- M[c1:c2, 1:n]-means
    #         #     }
    #         #     Pi <- c(Pi, 0.5 * means);gc()
    #         # }
    #         # print("Assignment DONE!")
    #         # if(is.null(SUM)){
    #         #     SUM <- sum(Pi * (1-Pi))
    #         # }
    #         # fl.suc <- flush(Z)
    #         # if(!fl.suc){ stop("flush failed\n") }
    #         # RR <- describe(Z); rm(list=c("Z", "Pi", "means")); gc()
    #         # Z <- attach.big.matrix(RR)
    #         # print("Computing Z'Z in big.matrix...")
    #         # K <- 0.5 * big.crossprod(Z)/SUM
    #         # rm(Z)
    #         # gc()
    #         # unlink(c(paste0("Z", memo, ".temp.bin"), paste0("Z", memo, ".temp.desc")), recursive = TRUE)
    #     }
    # )
    #print("K Preparation is Done!")
    logging.log("Deriving relationship matrix successfully", "\n", verbose = verbose); gc()
    return(K)
}#end of MVP.k.VanRaden function
