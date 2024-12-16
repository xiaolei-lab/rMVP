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

#ifndef MVP_OMP_H_
#define MVP_OMP_H_

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <Rcpp.h>

static int omp_setup(int threads);
static inline int omp_setup(int threads=0) {
    int t = 1;
#ifdef _OPENMP
    if (threads == 0) {
        t = omp_get_num_procs() - 1;
        t = t > 0 ? t : 1;
    } else {
        t = threads > 0 ? threads : 1;
    }
    omp_set_num_threads(t);
#else
#endif
    return t;
}

#endif
