#include "rMVP.h"

// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::mat GInv(const arma::mat A){
	
	arma::mat ginv;
	if(A.n_rows == 1){
		ginv = 1 / A;
	}else{
		arma::mat U;
		arma::vec s;
		arma::mat V;
		double tol = sqrt(datum::eps);
		
		svd(U,s,V,A);
		U = conv_to<mat>::from(conj(conv_to<cx_mat>::from(U)));
		arma::vec sMax(2); sMax.fill(0);
		sMax[1] = tol * s[0];
		arma::uvec Positive = find(s > sMax.max());
		arma::mat Up = U.cols(Positive);
		Up.each_row() %= 1/s(Positive).t();
		ginv = V.cols(Positive) * Up.t();
	}
	return ginv;
}

template <typename T>
NumericVector getRow(XPtr<BigMatrix> pMat, const int row){
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = pMat->ncol();

	NumericVector snp(ind);

	for(int i = 0; i < ind; i++){
		snp[i] = genomat[i][row];
	}

	return snp;
}

// [[Rcpp::export]]
NumericVector getRow(SEXP pBigMat, const int row){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return getRow<char>(xpMat, row);
	case 2:
		return getRow<short>(xpMat, row);
	case 4:
		return getRow<int>(xpMat, row);
	case 8:
		return getRow<double>(xpMat, row);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP glm_c(const arma::vec &y, const arma::mat &X, const arma::mat & iXX, XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	
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

	int q0 = X.n_cols;
	if(y.n_elem != n)	throw Rcpp::exception("number of individuals not match!");

	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

	// arma::mat iXX = GInv(X.t() * X);
	arma::mat xy = X.t() * y;
	double yy = sum(y % y);
	arma::mat res(m, 1 + 1 + 1 + q0);
	arma::mat iXXs(q0 + 1, q0 + 1);

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
			if(_marker_ind.is_empty()){
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[(i_marker + l)][k];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[k][(i_marker + l)];
						}
					}
				}
			}else{
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[_marker_ind[(i_marker + l)]][k];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[k][_marker_ind[(i_marker + l)]];
						}
					}
				}
			}
		}else{
			if(_marker_ind.is_empty()){
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[(i_marker + l)][_geno_ind[k]];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[_geno_ind[k]][(i_marker + l)];
						}
					}
				}
			}else{
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[_marker_ind[(i_marker + l)]][_geno_ind[k]];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[_geno_ind[k]][_marker_ind[(i_marker + l)]];
						}
					}
				}
			}
		}

		#pragma omp parallel for firstprivate(iXXs)
		for(int l = 0; l < cnt; l++){
			double sy = sum(Z_buffer.col(l) % y);
			double ss = sum(Z_buffer.col(l) % Z_buffer.col(l));
			arma::mat xs = X.t() * Z_buffer.col(l);
			arma::mat B21 = xs.t() * iXX;
			double t2 = as_scalar(B21 * xs);
			double B22 = (ss - t2);
			double invB22;
			int df;
			if(B22 < 1e-8){
				invB22 = 0;
				df = n - q0;
			}else{
				invB22 = 1 / B22;
				df = n - q0 - 1;
			}
			arma::mat NeginvB22B21 = -1 * invB22 * B21;
			iXXs(q0, q0) = invB22;
			iXXs.submat(0, 0, q0 - 1, q0 - 1) = iXX + invB22 * B21.t() * B21;
			iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
			iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

			// statistics
			arma::mat rhs(xy.n_rows + 1, 1);
			rhs.rows(0, xy.n_rows - 1) = xy;
			rhs(xy.n_rows, 0) = sy;
			arma::mat beta = iXXs * rhs;

			double ve = (yy - as_scalar(beta.t() * rhs)) / df;
			arma::vec se(q0 + 1);
			arma::vec pvalue(q0 + 1);
			for(int ff = 0; ff < (q0 + 1); ff++){
				se[ff] = sqrt(iXXs(ff, ff) * ve);
				pvalue[ff] = 2 * R::pt(abs(beta[ff] / se[ff]), df, false, false);
				res(i_marker + l, ff + 2) = pvalue[ff];
			}

			if(invB22 == 0){
				beta[q0] = NA_REAL;
				se[q0] = NA_REAL;
				res(i_marker + l, q0) = NA_REAL;
			}
			res(i_marker + l, 0) = beta[q0];
			res(i_marker + l, 1) = se[q0]; 
		}

		i = j;
		i_marker += cnt;
		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();
	return wrap(res);
}

// [[Rcpp::export]]
SEXP glm_c(const arma::vec & y, const arma::mat & X, const arma::mat & iXX, SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return glm_c<char>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 2:
		return glm_c<short>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 4:
		return glm_c<int>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 8:
		return glm_c<double>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP mlm_c(const arma::vec & y, const arma::mat & X, const arma::mat & U, const double vgs, XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	
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

	int q0 = X.n_cols;
	if(y.n_elem != n)	throw Rcpp::exception("number of individuals not match.!");

	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

	arma::mat Uy = U.t() * y;
	arma::mat UX = U.t() * X;
	arma::mat UXUy = UX.t() * Uy;
	arma::mat iUXUX = GInv(UX.t() * UX);
	
	arma::mat res(m, 3);		
	arma::mat iXXs(q0 + 1, q0 + 1);

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
			if(_marker_ind.is_empty()){
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[(i_marker + l)][k];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[k][(i_marker + l)];
						}
					}
				}
			}else{
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[_marker_ind[(i_marker + l)]][k];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[k][_marker_ind[(i_marker + l)]];
						}
					}
				}
			}
		}else{
			if(_marker_ind.is_empty()){
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[(i_marker + l)][_geno_ind[k]];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[_geno_ind[k]][(i_marker + l)];
						}
					}
				}
			}else{
				if(marker_bycol){
					#pragma omp parallel for
					for(int l = 0; l < cnt; l++){
						for(int k = 0; k < n; k++){
							Z_buffer(k, l) = (double)genomat[_marker_ind[(i_marker + l)]][_geno_ind[k]];
						}
					}
				}else{
					#pragma omp parallel for
					for(int k = 0; k < n; k++){
						for(int l = 0; l < cnt; l++){
							Z_buffer(k, l) = (double)genomat[_geno_ind[k]][_marker_ind[(i_marker + l)]];
						}
					}
				}
			}
		}

		#pragma omp parallel for firstprivate(iXXs)
		for(int l = 0; l < cnt; l++){
			arma::mat Us = U.t() * Z_buffer.col(l);
			arma::mat UXUs = UX.t() * Us;
			double UsUs = as_scalar(Us.t() * Us);
			double UsUy = as_scalar(Us.t() * Uy);
			double B22 = UsUs - as_scalar(UXUs.t() * iUXUX * UXUs);
			double invB22 = 1 / B22;
			arma::mat B21 = UXUs.t() * iUXUX;
			arma::mat NeginvB22B21 = -1 * invB22 * B21;
			
			iXXs(q0, q0)=invB22;
			iXXs.submat(0, 0, q0 - 1, q0 - 1) = iUXUX + invB22 * B21.t() * B21;
			iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
			iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

			// statistics
			arma::mat rhs(UXUy.n_rows + 1, 1);
			rhs.rows(0, UXUy.n_rows - 1) = UXUy;
			rhs(UXUy.n_rows, 0) = UsUy;
			arma::mat beta = iXXs * rhs;
			int df = n - q0 - 1;

			res(i_marker + l, 0) = beta(q0, 0);
			res(i_marker + l, 1) = sqrt(iXXs(q0, q0) * vgs); 
			res(i_marker + l, 2) = 2 * R::pt(abs(res(i_marker + l, 0) / res(i_marker + l, 1)), df, false, false);
		}

		i = j;
		i_marker += cnt;
		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();
	return wrap(res);
}

// [[Rcpp::export]]
SEXP mlm_c(const arma::vec & y, const arma::mat & X, const arma::mat & U, const double vgs, SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return mlm_c<char>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 2:
		return mlm_c<short>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 4:
		return mlm_c<int>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 8:
		return mlm_c<double>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
