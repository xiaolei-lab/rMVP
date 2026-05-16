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


# Internal FaSTLMM core — shared by MVP.FaSTLMM.LL and FarmCPU.FaSTLMM.LL.
# Returns list(beta, delta, LL, vg, ve).
.fastlmm_core <- function(pheno, snp.pool, X0=NULL, ncpus=2) {
    y <- pheno
    deltaExpStart <- .FASTLMM_DELTA_EXP_START
    deltaExpEnd   <- .FASTLMM_DELTA_EXP_END
    snp.pool <- snp.pool[,]
    if (!is.null(snp.pool) && any(apply(snp.pool, 2, var) == 0)) {
        deltaExpStart <- .FASTLMM_DELTA_EXP_DEGENERATE
        deltaExpEnd   <- deltaExpStart
    }
    if (is.null(X0)) X0 <- matrix(1, nrow(snp.pool), 1)
    X <- X0

    # SVD of snp.pool
    K.X.svd <- svd(snp.pool)
    d <- K.X.svd$d
    d <- d[d > .FASTLMM_SVD_THRESHOLD]
    d <- d^2
    U1 <- K.X.svd$u
    U1 <- U1[, 1:length(d)]
    if (is.null(dim(U1))) U1 <- matrix(U1, ncol=1)
    n <- nrow(U1)

    U1TX <- crossprod(U1, X)
    U1TY <- crossprod(U1, y)
    yU1TY <- y - U1 %*% U1TY
    XU1TX <- X - U1 %*% U1TX
    IU <- -tcrossprod(U1)
    diag(IU) <- rep(1, n) + diag(IU)
    IUX <- crossprod(IU, X)
    IUY <- crossprod(IU, y)

    delta.range <- seq(deltaExpStart, deltaExpEnd, by=.FASTLMM_DELTA_EXP_STEP)
    m <- length(delta.range)

    beta.optimize.parallel <- function(ii) {
        delta <- exp(delta.range[ii])
        dInv  <- 1 / (d + delta)

        beta1 <- crossprod(sweep(U1TX, 1, sqrt(dInv), "*"))
        beta2 <- crossprod(IUX) / delta
        beta3 <- crossprod(U1TX, U1TY * dInv)
        beta4 <- crossprod(IUX, IUY) / delta
        beta  <- crossprod(.safe_solve(beta1 + beta2), beta3 + beta4)

        part1 <- -0.5 * (n * log(.FASTLMM_2PI) + sum(log(d + delta)) + (n - length(d)) * log(delta))

        part221 <- sum((U1TY  - U1TX  %*% beta)^2 * dInv)
        part222 <- sum((yU1TY - XU1TX %*% beta)^2) / delta
        part2   <- -0.5 * (n + n * log((part221 + part222) / n))

        list(beta=beta, delta=delta, LL=part1 + part2)
    }

    llresults <- lapply(1:m, beta.optimize.parallel)

    LL_values <- sapply(llresults, function(x) x$LL)
    best_idx  <- which.max(LL_values)
    beta  <- llresults[[best_idx]]$beta
    delta <- llresults[[best_idx]]$delta
    LL    <- llresults[[best_idx]]$LL

    # Vectorized variance component estimates
    sigma_a1 <- sum((U1TY - U1TX %*% beta)^2 / (d + delta))
    sigma_a2 <- sum((IUY  - IUX  %*% beta)^2) / delta
    sigma_a  <- (sigma_a1 + sigma_a2) / n
    sigma_e  <- delta * sigma_a

    list(beta=beta, delta=delta, LL=LL, vg=sigma_a, ve=sigma_e)
}


#' Evaluation of the maximum likelihood using FaST-LMM method
#'
#' Last update: January 11, 2017
#'
#' @author Xiaolei Liu (modified)
#'
#' @param pheno a two-column phenotype matrix
#' @param snp.pool matrix for pseudo QTNs
#' @param X0 covariates matrix
#' @param ncpus number of threads used for parallel computation
#'
#' @return
#' Output: beta - beta effect
#' Output: delta - delta value
#' Output: LL - log-likelihood
#' Output: vg - genetic variance
#' Output: ve - residual variance
#'
`MVP.FaSTLMM.LL` <- function(pheno, snp.pool, X0=NULL, ncpus=2) {
    .fastlmm_core(pheno=pheno, snp.pool=snp.pool, X0=X0, ncpus=ncpus)
}
