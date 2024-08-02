#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include "rMVP.h"
#include <R_ext/Print.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

template <typename T>
arma::vec BigRowMean(XPtr<BigMatrix> pMat, int threads = 0){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->ncol();
	int j, k, m = pMat->nrow();
	double p1 = 0.0;
	arma::vec mean(m);

	#pragma omp parallel for private(p1, k)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			p1 += bigm[k][j];
		}
		mean[j] = p1 / ind;
	}

	return mean;
}


arma::vec BigRowMean(SEXP pBigMat, int threads = 0){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return BigRowMean<char>(xpMat, threads);
	case 2:
		return BigRowMean<short>(xpMat, threads);
	case 4:
		return BigRowMean<int>(xpMat, threads);
	case 8:
		return BigRowMean<double>(xpMat, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}


template <typename T>
SEXP kin_cal_m(XPtr<BigMatrix> pMat, int threads = 0, bool verbose = true){

    omp_setup(threads);

	if(verbose)
		Rcout << "Computing GRM under mode: Memory" << endl;

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = pMat->ncol();
	int m = pMat->nrow();
	int i = 0, j = 0, k = 0;
	MinimalProgressBar pb;

	arma::vec Mean = BigRowMean(pMat, threads);
	double SUM = sum((0.5 * Mean) % (1 - 0.5 * Mean));

	arma::mat kin(n, n);
	arma::vec coli(m);
	arma::vec colj(m);

	Progress p(n, verbose, pb);

	if(verbose)
		Rcout << "Scale the genotype matrix and compute Z'Z" << endl;

	#pragma omp parallel for schedule(dynamic) firstprivate(coli, colj) private(i, j, k) 
	for(i = 0; i < n; i++){
		for(k = 0; k < m; k++){
			coli[k] = bigm[i][k] - Mean[k];
		}
		if ( ! Progress::check_abort() ) {
			p.increment();
			for(j = i; j < n; j++){
				for(k = 0; k < m; k++){
					colj[k] = bigm[j][k] - Mean[k];
				}
				kin(i, j) = kin(j, i) = 0.5 * sum(coli % colj) / SUM;
			}
		}
	}

	return Rcpp::wrap(kin);
}


// [[Rcpp::export]]
SEXP kin_cal_m(SEXP pBigMat, int threads = 0, bool verbose = true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return kin_cal_m<char>(xpMat, threads, verbose);
	case 2:
		return kin_cal_m<short>(xpMat, threads, verbose);
	case 4:
		return kin_cal_m<int>(xpMat, threads, verbose);
	case 8:
		return kin_cal_m<double>(xpMat, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}


template <typename T>
SEXP kin_cal_s(XPtr<BigMatrix> pMat, int threads = 0, bool mkl = false, bool verbose = true){

    omp_setup(threads);

	if(verbose)
		Rcout << "Computing GRM under mode: Speed" << endl;

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	#ifdef _OPENMP
	#else
		if(!mkl)
			mkl = true;
	#endif
	if(threads == 1)
		mkl = true;

	int n = pMat->ncol();
	int m = pMat->nrow();
	int i = 0, j = 0;
	MinimalProgressBar pb;

	arma::vec Mean = BigRowMean(pMat, threads);
	double SUM = sum((0.5 * Mean) % (1 - 0.5 * Mean));

	arma::mat kin(n, n);
	arma::mat geno(m, n);

	if(verbose)
		Rcout << "Scale the genotype matrix" << endl;

	#pragma omp parallel for schedule(dynamic) private(i, j)
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			geno(j, i) = bigm[i][j] - Mean[j];
		}
	}

	if(verbose)
		Rcout << "Computing Z'Z" << endl;

	if(mkl){
		kin = geno.t() * geno / SUM / 2;
	}else{

		Progress p(n, verbose, pb);
		arma::colvec coli;

		#pragma omp parallel for schedule(dynamic) private(i, j, coli)
		for(i = 0; i < n; i++){
			coli = geno.col(i);
			if ( ! Progress::check_abort() ) {
				p.increment();
				for(j = i; j < n; j++){
					kin(j, i) = kin(i, j) = 0.5 * sum(coli % geno.col(j)) / SUM;
				}
			}
		}

	}

	return Rcpp::wrap(kin);
}


// [[Rcpp::export]]
SEXP kin_cal_s(SEXP pBigMat, int threads = 0, bool mkl = false, bool verbose = true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return kin_cal_s<char>(xpMat, threads, mkl, verbose);
	case 2:
		return kin_cal_s<short>(xpMat, threads, mkl, verbose);
	case 4:
		return kin_cal_s<int>(xpMat, threads, mkl, verbose);
	case 8:
		return kin_cal_s<double>(xpMat, threads, mkl, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}


template <typename T>
SEXP kin_cal(XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, int threads = 0, size_t step = 5000, bool mkl = false, bool verbose = true){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	#ifdef _OPENMP
	#else
		if(!mkl)
			mkl = true;
	#endif
	if(threads == 1)
		mkl = true;

	uvec _geno_ind;
	if(geno_ind.isNotNull()){
        _geno_ind = as<uvec>(geno_ind) - 1;
    }else{
		_geno_ind = regspace<uvec>(0, pMat->ncol() - 1);
	}

	int n = _geno_ind.n_elem;
	int m = pMat->nrow();

	arma::vec means(m);
	#pragma omp parallel for
	for (size_t j = 0; j < m; j++){
		double p1 = 0.0;
		for(size_t k = 0; k < n; k++){
			p1 += bigm[_geno_ind[k]][j];
		}
		means[j] = p1 / n;
	}

	arma::mat kin = zeros<mat>(n, n);
	arma::mat Z_buffer(step, n, fill::none);

	size_t i = 0, j = 0;
	size_t i_marker = 0;
	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);
	for (;i < m;) {
		
		int cnt = 0;
		for (; j < m && cnt < step; j++)
		{
			cnt++;
		}

		if (cnt != step) {
			Z_buffer.set_size(cnt, n);
		}

		#pragma omp parallel for
		for(size_t k = 0; k < n; k++){
			for(size_t l = 0; l < cnt; l++){
				Z_buffer(l, k) = bigm[_geno_ind[k]][i_marker + l] - means[i_marker + l];
			}
		}

		if(mkl){
			double alp = 1.0;
			double beta = 1.0;
			char uplo = 'L';
			dsyrk_(&uplo, "T", &n, &cnt, &alp, Z_buffer.memptr(), &cnt, &beta, kin.memptr(), &n);
		}else{
			arma::colvec coli;
			#pragma omp parallel for schedule(dynamic) private(coli)
			for(size_t k = 0; k < n; k++){
				coli = Z_buffer.col(k);
				for(size_t l = k; l < n; l++){
					kin(l, k) += sum(coli % Z_buffer.col(l));
				}
			}
		}
		i = j;
		i_marker += cnt;

		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();

	#pragma omp parallel for schedule(dynamic)
	for (uword j = 0; j < kin.n_cols; j++) {
		for (uword i = (j + 1); i < kin.n_cols; i++) {
			kin(j, i) = kin(i, j);
		}
	}
	kin /= arma::mean(kin.diag());

	return Rcpp::wrap(kin);
}

// [[Rcpp::export]]
SEXP kin_cal(SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, int threads = 0, size_t step = 10000, bool mkl = false, bool verbose = true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return kin_cal<char>(xpMat, geno_ind, threads, step, mkl, verbose);
	case 2:
		return kin_cal<short>(xpMat, geno_ind, threads, step, mkl, verbose);
	case 4:
		return kin_cal<int>(xpMat, geno_ind, threads, step, mkl, verbose);
	case 8:
		return kin_cal<double>(xpMat, geno_ind, threads, step, mkl, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
