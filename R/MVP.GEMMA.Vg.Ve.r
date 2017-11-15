MVP.GEMMA.Vg.Ve <- function(y, X, K, rtol=1e-6, atol=1e-8, ctol=1e-8, root=FALSE){
##########################################################################################################
#Object: To estimate variance component using HE regression
#Input:
#y: phenotype
#X: genotype
#K: kinship matrix
#rtol, atol, ctol, root: parameters for HE regression, no changes is recommended
#Output: vg, ve, and delta
#Author: Xiang Zhou
#Translated from GEMMA (C++) to R by: Haohao Zhang
#Modified by: Lilin Yin and Xiaolei Liu
#Build date: Feb 2, 2017
#Last update: Feb 2, 2017
##########################################################################################################

	#try(setMKLthreads(1),silent = TRUE)
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

    # multiroot
	if(root){
		mult.res = multiroot(LogRL_dev1, start = log_sigma2, parms=list(y = y,K = K,X = X),rtol = rtol, atol = atol, ctol = ctol, maxiter=1000)
		vg = exp(mult.res$root[1])
		ve = exp(mult.res$root[2])
		delta = ve / vg
	}else{
		vg = exp(log_sigma2[1])
		ve = exp(log_sigma2[2])
		delta = ve / vg
	}
	print("Variance Components Estimation is Done!")
	return(list(vg=vg, ve=ve, delta=delta))
}
