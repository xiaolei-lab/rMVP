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


#' To estimate variance component using HE regression
#' 
#' Build date: Feb 2, 2017
#' Last update: Feb 2, 2019
#' 
#' @author Translated from C++(GEMMA, Xiang Zhou) to R by: Haohao Zhang
#' 
#' @param y phenotype
#' @param X genotype
#' @param K kinship matrix

#' @return vg, ve, and delta
#' @export
#'
#' @examples
#' \donttest{
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' 
#' K <- MVP.K.VanRaden(genotype, cpu=1)
#' vc <- MVP.HE.Vg.Ve(y=phenotype[,2], X=matrix(1, nrow(phenotype)), K=K)
#' print(vc)
#' }
#' 
MVP.HE.Vg.Ve <- function(y, X, K) {
    # NA in phenotype
    idx <- !is.na(y)
    y <- y[idx]
    X <- as.matrix(X[idx, ])
    K <- K[idx, idx]
    
    K = K[]
    n = nrow(K)

    #center the Kinship
    w = rep(1, n)
    Gw = rowSums(K)
    alpha = -1.0 / nrow(K)
    #K = K + alpha * (w %*% t(Gw) + Gw %*% t(w))
    K = K + alpha * (tcrossprod(w, Gw) + tcrossprod(Gw, w))
    
    CalcVChe <- function(y, X, K){
        n = nrow(K)
        r = n / (n - ncol(X))
        v_traceG = mean(diag(K))
        
        # center and scale K by X
        WtW = crossprod(X)
        WtWi = solve(WtW)
        WtWiWt = tcrossprod(WtWi, X)
        GW = K %*% X
        Gtmp = GW %*% WtWiWt
        K = K - Gtmp
        Gtmp = t(Gtmp)
        K = K - Gtmp
        WtGW = crossprod(X, GW)
        GW = crossprod(WtWiWt, WtGW)
        Gtmp = GW %*% WtWiWt
        K = K + Gtmp
        d = mean(diag(K))
        traceG_new = d
        if (d != 0) {
            K = K * 1 / d
        }

        # center y by X, and standardize it to have variance 1 (t(y)%*%y/n=1)
        Wty = crossprod(X, y)
        WtWiWty = solve(WtW, Wty)
        y_scale = y - X %*% WtWiWty
        
        VectorVar <- function(x) {
            m = mean(x)
            m2 = sum(x^2) / length(x)
            return (m2 - m * m)
        }

        StandardizeVector <- function(x) {
            m = mean(x)
            v = sum((x - m)^2) / length(x)
            v = v - m * m
            x = (x - m) / sqrt(v)
            return (x)
        }
        
        var_y = VectorVar(y)
        var_y_new = VectorVar(y_scale)

        y_scale = StandardizeVector(y_scale)

        # compute Kry, which is used for confidence interval; also compute q_vec (*n^2)
        Kry = K %*% y_scale - r * y_scale
        q_vec = crossprod(Kry, y_scale)

        # compuate yKrKKry, which is used later for confidence interval
        KKry = K %*% Kry
        yKrKKry = crossprod(Kry, KKry)
        d = crossprod(Kry)
        yKrKKry = c(yKrKKry, d)

        # compute Sij (*n^2)
        tr = sum(K^2)
        S_mat = tr - r * n

        # compute S^{-1}q
        Si_mat = 1/S_mat

        # compute pve (on the transformed scale)
        pve = Si_mat * q_vec

        # compute q_var (*n^4)
        s = 1
        qvar_mat = yKrKKry[1] * pve
        s = s - pve
        tmp_mat = yKrKKry[2] * s
        qvar_mat = (qvar_mat + tmp_mat) * 2.0

        # compute S^{-1}var_qS^{-1}
        Var_mat = Si_mat * qvar_mat * Si_mat

        # transform pve back to the original scale and save data
        s = 1
        vyNtgN = var_y_new / traceG_new
        vtgvy = v_traceG / var_y
        v_sigma2 = pve * vyNtgN
        v_pve = pve * (vyNtgN) * (vtgvy)
        s = s - pve
        pve_total = pve * (vyNtgN) * (vtgvy)

        d = sqrt(Var_mat)
        v_se_sigma2 = d * vyNtgN
        v_se_pve = d * (vyNtgN) * (vtgvy)
        se_pve_total = Var_mat * (vyNtgN) * (vtgvy) * (vyNtgN) * (vtgvy)
        
        v_sigma2 = c(v_sigma2, s * r * var_y_new)
        v_se_sigma2 = c(v_se_sigma2, sqrt(Var_mat) * r * var_y_new)
        se_pve_total = sqrt(se_pve_total)

        return(v_sigma2)
    }

    # initialize sigma2/log_sigma2
    v_sigma2 = CalcVChe(y=y, X=X, K=K)
    
    log_sigma2 = NULL
    for (i in 1:length(v_sigma2)) {
        if (v_sigma2[i] <= 0){
            log_sigma2[i] = log(0.1)
        } else {
            log_sigma2[i] = log(v_sigma2[i])
        }
    }
    
    LogRL_dev1 <- function(log_sigma2, parms) {
        
        y = parms$y
        X = parms$X
        K = parms$K
        n = nrow(K)
        
        P = K * exp(log_sigma2[1]) + diag(n) * exp(log_sigma2[2])

        # calculate H^{-1}
        P = solve(P)

        # calculate P=H^{-1}-H^{-1}X(X^TH^{-1}X)^{-1}X^TH^{-1}
        HiW = P %*% X
        WtHiW = crossprod(X, HiW)
        WtHiWi = solve(WtHiW)

        WtHiWiWtHi = tcrossprod(WtHiWi, HiW)
        P = P - HiW %*% WtHiWiWtHi

        # calculate Py, KPy, PKPy
        Py = P %*% matrix(y)
        KPy = K %*% Py
        
        # calculate dev1=-0.5*trace(PK_i)+0.5*yPKPy
        c(para1 = (-0.5 * sum(t(P) * K) + 0.5 * crossprod(Py, KPy)) * exp(log_sigma2[1]),
            para2 = (-0.5 * sum(diag(P)) + 0.5 * crossprod(Py)) * exp(log_sigma2[2]))
    }

    vg = exp(log_sigma2[1])
    ve = exp(log_sigma2[2])
    delta = ve / vg
    return(list(vg=vg, ve=ve, delta=delta))
}
