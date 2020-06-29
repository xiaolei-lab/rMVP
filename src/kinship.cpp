#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include "mvp_omp.h"
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


class MinimalProgressBar: public ProgressBar{
	public:
	MinimalProgressBar()  {
		_finalized = false;
	}
	~MinimalProgressBar() {}
	void display() {}
	void update(float progress) {
		if (_finalized) return;
		REprintf("\r");
		REprintf("Calculating in process...(finished %.2f%)", progress * 100);
	}
	void end_display() {
	if (_finalized) return;
		REprintf("\r");
		
		REprintf("Calculating in process...(finished 100.00%)");
		REprintf("\n");
		_finalized = true;
	}
	private:
	bool _finalized;
};


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
	double p12 = 0.0;
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
