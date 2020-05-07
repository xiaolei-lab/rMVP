// Copyright (C) 2016-2019 by Xiaolei Team
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include "mvp_omp.h"
#include <progress.hpp>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::depends(RcppProgress)]]

template <typename T>
void impute_marker(XPtr<BigMatrix> pMat, int threads=0, bool verbose=true) {
    omp_setup(threads);
    Progress progress(pMat->nrow(), verbose);
    
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    const size_t n = pMat->ncol();
    const size_t m = pMat->nrow();
    
    // loop marker
    #pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        std::vector<size_t> na_index = {};;
        size_t counts[3] = { 0 };
        
        // count allele, record missing index 
        for (size_t j = 0; j < n; j++) {
            switch(int(mat[j][i])) {
            case 0: counts[0]++; break;
            case 1: counts[1]++; break;
            case 2: counts[2]++; break;
            default: na_index.push_back(j);
            }
        }
        
        // find major allele
        T major = counts[2] > counts[1] ? (counts[2] > counts[0] ? 2 : 0) : (counts[1] > counts[0] ? 1 : 0);
        
        // impute
        for (auto&& x: na_index) {
            mat[x][i] = major;   
        }
        progress.increment();
    }
}

// [[Rcpp::export]]
void impute_marker(SEXP pBigMat, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return impute_marker<char>(xpMat, threads, verbose);
    case 2:
        return impute_marker<short>(xpMat, threads, verbose);
    case 4:
        return impute_marker<int>(xpMat, threads, verbose);
    case 8:
        return impute_marker<double>(xpMat, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

template <typename T>
bool hasNA(XPtr<BigMatrix> pMat, double NA_C) {
    size_t m = pMat->nrow();
    size_t n = pMat->ncol();

    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < m; i++) {
            if (mat[j][i] == NA_C) {
                return true;
            }
        }
    }
    return false;
}

// [[Rcpp::export]]
bool hasNA(SEXP pBigMat) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return hasNA<char>(xpMat, NA_CHAR);
    case 2:
        return hasNA<short>(xpMat, NA_SHORT);
    case 4:
        return hasNA<int>(xpMat, NA_INTEGER);
    case 8:
        return hasNA<double>(xpMat, NA_REAL);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

/*** R

*/
