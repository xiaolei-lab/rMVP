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
arma::vec BigRowMean(XPtr<BigMatrix> pMat, bool marker_bycol = true, size_t step = 10000, int threads = 0, const Nullable<arma::uvec> geno_ind = R_NilValue, const bool verbose = true){

    omp_setup(threads);
	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n;
	int m = marker_bycol ? pMat->ncol() : pMat->nrow();
	arma::vec mean(m, fill::zeros);

	uvec _geno_ind;
	if(geno_ind.isNotNull()){
        _geno_ind = as<uvec>(geno_ind) - 1;
        n = _geno_ind.n_elem;
    }else{
		n = marker_bycol ? pMat->nrow() : pMat->ncol();
	}

	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

	arma::mat Z_buffer(n, step, fill::none);
	int i = 0, j = 0;
	int i_marker = 0;
	for (;i < m;) {
		
		int cnt = 0;
		for (; j < m && cnt < step; j++)
		{
			cnt++;
		}

		if (cnt != step) {
			Z_buffer.set_size(n, cnt);
		}

		if(_geno_ind.is_empty()){
			if(marker_bycol){
				#pragma omp parallel for
				for(int l = 0; l < cnt; l++){
					for(int k = 0; k < n; k++){
						T elem = bigm[(i_marker + l)][k];
						Z_buffer(k, l) = (isna(elem) ? datum::nan : (double)elem);
					}
				}
			}else{
				#pragma omp parallel for
				for(int k = 0; k < n; k++){
					for(int l = 0; l < cnt; l++){
						T elem = bigm[k][(i_marker + l)];
						Z_buffer(k, l) = (isna(elem) ? datum::nan : (double)elem);
					}
				}
			}
		}else{
			if(marker_bycol){
				#pragma omp parallel for
				for(int l = 0; l < cnt; l++){
					for(int k = 0; k < n; k++){
						T elem = bigm[(i_marker + l)][_geno_ind[k]];
						Z_buffer(k, l) = (isna(elem) ? datum::nan : (double)elem);
					}
				}
			}else{
				#pragma omp parallel for
				for(int k = 0; k < n; k++){
					for(int l = 0; l < cnt; l++){
						T elem = bigm[_geno_ind[k]][(i_marker + l)];
						Z_buffer(k, l) = (isna(elem) ? datum::nan : (double)elem);
					}
				}
			}
		}

		#pragma omp parallel for
		for(int l = 0; l < cnt; l++){
			mean[i_marker + l] = arma::mean(Z_buffer.col(l));
		}
		if(mean.subvec(i_marker, i_marker + cnt - 1).has_nan())	throw Rcpp::exception("NA is not allowed in genotype, use 'MVP.Data.impute' to impute!");

		i = j;
		i_marker += cnt;
		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();

	return mean;
}

// [[Rcpp::export]]
arma::vec BigRowMean(SEXP pBigMat, bool marker_bycol = true, size_t step = 10000, int threads = 0, const Nullable<arma::uvec> geno_ind = R_NilValue, const bool verbose = true){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	DISPATCH_MATRIX_TYPE(BigRowMean, xpMat, marker_bycol, step, threads, geno_ind, verbose);
}


template <typename T>
SEXP kin_cal(XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const Nullable<arma::vec> marker_freq = R_NilValue, const bool marker_bycol = true, int threads = 0, int step = 5000, bool mkl = false, bool verbose = true){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	#ifdef _OPENMP
	#else
		if(!mkl)
			mkl = true;
	#endif
	if(threads == 1)
		mkl = true;

	int n;
	uvec _geno_ind;
	if(geno_ind.isNotNull()){
        _geno_ind = as<uvec>(geno_ind) - 1;
		n = _geno_ind.n_elem;
    }else{
		n = marker_bycol ? pMat->nrow() : pMat->ncol();
	}

	int m;
	uvec _marker_ind;
	if(marker_ind.isNotNull()){
        _marker_ind = as<uvec>(marker_ind) - 1;
		m = _marker_ind.n_elem;
    }else{
		m = marker_bycol ? pMat->ncol() : pMat->nrow();
	}

	vec means;
	if(marker_freq.isNotNull()){
		means = as<vec>(marker_freq) * 2;
		if(means.has_nan())	throw Rcpp::exception("NA is not allowed in allele frequency!");
	}

	arma::mat kin = zeros<mat>(n, n);
	arma::mat Z_buffer(step, n, fill::none);

	int i = 0, j = 0;
	int i_marker = 0;
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

		fill_geno_buffer(Z_buffer, bigm, cnt, i_marker, _geno_ind, _marker_ind, marker_bycol, false);

		if(marker_freq.isNotNull()){
			Z_buffer.each_col() -= means.subvec(i_marker, i_marker + cnt - 1);
		}else{
			means = mean(Z_buffer, 1);
			Z_buffer.each_col() -= means;
		}

		if(mkl){
			double alp = 1.0;
			double beta = 1.0;
			char uplo = 'L';
			dsyrk_(&uplo, "T", &n, &cnt, &alp, Z_buffer.memptr(), &cnt, &beta, kin.memptr(), &n);
		}else{
			arma::colvec coli;
			#pragma omp parallel for schedule(dynamic) private(coli)
			for(int k = 0; k < n; k++){
				coli = Z_buffer.col(k);
				for(int l = k; l < n; l++){
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
SEXP kin_cal(SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const Nullable<arma::vec> marker_freq = R_NilValue, const bool marker_bycol = true, int threads = 0, int step = 10000, bool mkl = false, bool verbose = true){

	XPtr<BigMatrix> xpMat(pBigMat);

	DISPATCH_MATRIX_TYPE(kin_cal, xpMat, geno_ind, marker_ind, marker_freq, marker_bycol, threads, step, mkl, verbose);
}
