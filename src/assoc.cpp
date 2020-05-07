#include <RcppArmadillo.h>
#include "mvp_omp.h"
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

class MinimalProgressBar_plus: public ProgressBar{
	public:
	MinimalProgressBar_plus()  {
		_finalized = false;
	}

	~MinimalProgressBar_plus() {}

	std::string _time_to_string(double seconds) {
	  
	    int time = (int) seconds;
	  
	    int hour = 0;
	    int min = 0;
	    int sec = 0;
	  
	    hour = time / 3600;
	    time = time % 3600;
	    min = time / 60;
	    time = time % 60;
	    sec = time;
	  
	    std::stringstream time_strs;
	    time_strs << "TimeLeft: ";
	    if (hour != 0) time_strs << hour << "h";
	    if (hour != 0 || min != 0) time_strs << min << "m";
	    time_strs << sec << "s";
	    std::string time_str = time_strs.str();
	  
	    return time_str;
	}

	int _compute_nb_ticks(float progress) {
	    return int(progress * _max_ticks);
	}

	std::string _construct_ticks_display_string(int nb) {
	    std::stringstream ticks_strs;
	    for (int i = 1; i <= _max_ticks; ++i) {
	    	if(i < 4){
	    		ticks_strs << ">";
	    	} else if (i < nb) {
	            ticks_strs << "-";
	        } else if(i == nb) {
	        	ticks_strs << ">";
	        } else {
	            ticks_strs << " ";
	        }
	    }
	    std::string tick_space_string = ticks_strs.str();
	    return tick_space_string;
	}

	void end_display() {
      update(1);
    }

    void _finalize_display() {
      if (_finalized) return;
      REprintf("\n");
      _finalized = true;
    }

	void display() {
		REprintf("\r");
	}

	void update(float progress) {
	  
	    // stop if already finalized
	    if (_finalized) return;
	  
	    // start time measurement when update() is called the first time
	    if (_timer_flag) {
	        _timer_flag = false;
	        // measure start time
	        time(&start);
	    } else {
	    
	    	int nb_ticks = _compute_nb_ticks(progress);
		    int delta = nb_ticks - _ticks_displayed;
		    if (delta > 0) {
		    	_ticks_displayed = nb_ticks;
		        std::string cur_display = _construct_ticks_display_string(nb_ticks);
		        
		        // measure current time
		        time(&end);
	    
		        // calculate passed time and remaining time (in seconds)
		        double pas_time = std::difftime(end, start);
		        double rem_time = (pas_time / progress) * (1 - progress);
		        if(rem_time < 1 && rem_time > 0.5)	rem_time = 1;

		    	// convert seconds to time string
		        std::string time_string = _time_to_string(rem_time);
		        
		        // ensure overwriting of old time info
		        int empty_length = time_string.length();
		        
	        	std::string empty_space;

		        // merge progress bar and time string
		        std::stringstream strs;
		        if(empty_length_p && abs(empty_length - empty_length_p)){
		        	empty_space = std::string(abs(empty_length - empty_length_p), ' ');
		        	strs << "[" << cur_display << "] " << time_string << empty_space;
		        }else{
		        	strs << "[" << cur_display << "] " << time_string;
		        }
		    	empty_length_p = empty_length;
		    	// strs << "[" << cur_display << "]";

		        std::string temp_str = strs.str();
		        char const* char_type = temp_str.c_str();
		    
		        // print: remove old display and replace with new
		        REprintf("\r");
		        REprintf("%s", char_type);

		    }
		    if (_ticks_displayed >= _max_ticks)
       			// end_display();
       			_finalize_display();
	    }
	}
	private: 
	int	empty_length_p = 0;
    int _max_ticks = 45;
    bool _finalized;
    bool _timer_flag = true;
    time_t start, end;
    int _ticks_displayed = 0;
};

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
SEXP glm_c(const arma::vec y, const arma::mat X, const arma::mat iXX, XPtr<BigMatrix> pMat, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = pMat->ncol();
	int mkr = pMat->nrow();
	int q0 = X.n_cols;

	int y_len = y.n_elem;
	if(y_len != ind)
		throw Rcpp::exception("number of individuals not match.!");

	MinimalProgressBar_plus pb;
	Progress progress(mkr, verbose, pb);

	// arma::mat iXX = GInv(X.t() * X);
	arma::mat xy = X.t() * y;
	double yy = dot(y, y);
	arma::mat res(mkr, 1 + 1 + q0);
	arma::vec snp(ind);
	arma::mat iXXs(q0 + 1, q0 + 1);

	#pragma omp parallel for schedule(dynamic) firstprivate(snp, iXXs)
	for(int i = 0; i < mkr; i++){

		for(int ii = 0; ii < ind; ii++){
			snp[ii] = genomat[ii][i];
		}
		
		double sy = dot(snp, y);
		double ss = dot(snp, snp);
		arma::mat xs = X.t() * snp;
		arma::mat B21 = xs.t() * iXX;
		double t2 = as_scalar(B21 * xs);
		double invB22 = 1 / (ss - t2);
		arma::mat NeginvB22B21 = -1 * invB22 * B21;

		iXXs(q0, q0)=invB22;

		iXXs.submat(0, 0, q0 - 1, q0 - 1) = iXX + invB22 * B21.t() * B21;
		iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
		iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

        // statistics
        arma::mat rhs(xy.n_rows + 1, 1);
        rhs.rows(0, xy.n_rows - 1) = xy;
        rhs(xy.n_rows, 0) = sy;
		arma::mat beta = iXXs * rhs;
        int df = ind - q0 - 1;
        double ve = (yy - as_scalar(beta.t() * rhs)) / df;
       	arma::vec se(q0 + 1);
        arma::vec pvalue(q0 + 1);
        for(int ff = 0; ff < (q0 + 1); ff++){
        	se[ff] = sqrt(iXXs(ff, ff) * ve);
        	pvalue[ff] = 2 * R::pt(abs(beta[ff] / se[ff]), df, false, false);
        	if(ff > 0)	res(i, ff + 1) = pvalue[ff];
        }
        res(i, 0) = beta[q0];
        res(i, 1) = se[q0]; 
        progress.increment();
	}

	return wrap(res);
}

// [[Rcpp::export]]
SEXP glm_c(const arma::vec y, const arma::mat X, const arma::mat iXX, SEXP pBigMat, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return glm_c<char>(y, X, iXX, xpMat, verbose, threads);
	case 2:
		return glm_c<short>(y, X, iXX, xpMat, verbose, threads);
	case 4:
		return glm_c<int>(y, X, iXX, xpMat, verbose, threads);
	case 8:
		return glm_c<double>(y, X, iXX, xpMat, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP mlm_c(const arma::vec y, const arma::mat X, const arma::mat U, const double vgs, XPtr<BigMatrix> pMat, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = pMat->ncol();
	int mkr = pMat->nrow();
	int q0 = X.n_cols;

	int y_len = y.n_elem;
	if(y_len != ind)
		throw Rcpp::exception("number of individuals not match.!");

	MinimalProgressBar_plus pb;
	Progress progress(mkr, verbose, pb);

	arma::mat Uy = U.t() * y;
	arma::mat UX = U.t() * X;
	arma::mat UXUy = UX.t() * Uy;
	arma::mat iUXUX = GInv(UX.t() * UX);
	
	arma::mat res(mkr, 3);		
	arma::vec snp(ind);
	arma::mat iXXs(q0 + 1, q0 + 1);

	#pragma omp parallel for schedule(dynamic) firstprivate(snp, iXXs)
	for(int i = 0; i < mkr; i++){

		for(int ii = 0; ii < ind; ii++){
			snp[ii] = genomat[ii][i];
		}
		arma::mat Us = U.t() * snp;
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
		int df = ind - q0 - 1;

        res(i, 0) = beta(q0, 0);
        res(i, 1) = sqrt(iXXs(q0, q0) * vgs); 
        res(i, 2) = 2 * R::pt(abs(res(i, 0) / res(i, 1)), df, false, false);
        progress.increment();
	}

	return wrap(res);
}

// [[Rcpp::export]]
SEXP mlm_c(const arma::vec y, const arma::mat X, const arma::mat U, const double vgs, SEXP pBigMat, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return mlm_c<char>(y, X, U, vgs, xpMat, verbose, threads);
	case 2:
		return mlm_c<short>(y, X, U, vgs, xpMat, verbose, threads);
	case 4:
		return mlm_c<int>(y, X, U, vgs, xpMat, verbose, threads);
	case 8:
		return mlm_c<double>(y, X, U, vgs, xpMat, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
