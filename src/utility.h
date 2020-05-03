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

#include <omp.h>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

static int omp_setup(int threads, bool verbose);
static inline int omp_setup(int threads=0, bool verbose=true) {
    int t = 1;
#ifdef _OPENMP
    if (threads == 0) {
        t = omp_get_num_procs() - 1;
        t = t > 0 ? t : 1;
    } else {
        t = threads > 0 ? threads : 1;
    }
    omp_set_num_threads(t);
    
    if (verbose)
        Rcpp::Rcerr << "Number of threads: " << omp_get_max_threads() << std::endl;
#else
    if (verbose)
        Rcpp::Rcerr << "Number of threads: 1 (No OpenMP detected)" << std::endl;
#endif
    return t;
}
