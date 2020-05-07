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
SEXP kin_cal(XPtr<BigMatrix> pMat, int threads = 0){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = pMat->ncol();
	int m = pMat->nrow();
	int i = 0, j = 0, k = 0;
	double p12 = 0.0;
	MinimalProgressBar pb;

	arma::vec Mean = BigRowMean(pMat, threads);
	double SUM = arma::dot((0.5 * Mean), (1 - 0.5 * Mean));

	arma::mat kin(n, n);

	Progress p(n, true, pb);

	#pragma omp parallel for schedule(dynamic) private(j, k, p12)
	for(i = 0; i < n; i++){
		if ( ! Progress::check_abort() ) {
			p.increment();
			for(j = i; j < n; j++){
				p12 = 0.0;
				for(k = 0; k < m; k++){
					p12 += (bigm[i][k] - Mean[k]) * (bigm[j][k] - Mean[k]);
				}
				kin(i, j) = 0.5 * p12 / SUM;
				kin(j, i) = kin(i, j);
			}
		}
	}

	// #pragma omp parallel for schedule(dynamic) private(j, p12)
	// std::vector<double> Means = as<std::vector<double> >(wrap(Mean));
	// for(i = 0; i < n; i++){
		// vector<int> m_i(bigm[i], bigm[i] + m);
		// for(j = i; j < n; j++){
			// vector<int> m_j(bigm[j], bigm[j] + m);
			// p12 = std::inner_product(m_i.begin(), m_i.end(), m_j.begin(), 0);
			// kin(i, j) = p12 / SUM;
			// kin(j, i) = kin(i, j);
		// }
	// }
	
	return Rcpp::wrap(kin);
}


// [[Rcpp::export]]
SEXP kin_cal(SEXP pBigMat, int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return kin_cal<char>(xpMat, threads);
	case 2:
		return kin_cal<short>(xpMat, threads);
	case 4:
		return kin_cal<int>(xpMat, threads);
	case 8:
		return kin_cal<double>(xpMat, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
