#ifndef RMVP_H_
#define RMVP_H_

#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS 1

#include <RcppArmadillo.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"
#include "mvp_omp.h"

// ---------------------------------------------------------------------------
// MinimalProgressBarBase — shared implementation for timed progress bars.
// Concrete subclasses must implement display() and _construct_ticks_display_string().
// ---------------------------------------------------------------------------
class MinimalProgressBarBase : public ProgressBar {
public:
    explicit MinimalProgressBarBase(int max_ticks)
        : _max_ticks(max_ticks), empty_length_p(0), _finalized(false),
          _timer_flag(true), _ticks_displayed(0) {}

    virtual ~MinimalProgressBarBase() {}

    virtual void display() = 0;

    void end_display() { update(1); }

    void update(float progress) {
        if (_finalized) return;

        if (_timer_flag) {
            _timer_flag = false;
            time(&start);
        } else {
            int nb_ticks = _compute_nb_ticks(progress);
            int delta = nb_ticks - _ticks_displayed;
            if (delta > 0) {
                _ticks_displayed = nb_ticks;
                std::string cur_display = _construct_ticks_display_string(nb_ticks);

                time(&end);
                double pas_time = std::difftime(end, start);
                double rem_time = (progress < 1.0 ? (pas_time / progress) * (1 - progress) : pas_time);
                if (rem_time < 1 && rem_time > 0.5) rem_time = 1;

                std::string time_string = _time_to_string(rem_time, progress);
                int empty_length = time_string.length();
                std::string empty_space;

                std::stringstream strs;
                if (empty_length_p && abs(empty_length - empty_length_p)) {
                    empty_space = std::string(abs(empty_length - empty_length_p), ' ');
                    strs << "[" << cur_display << "] " << time_string << empty_space;
                } else {
                    strs << "[" << cur_display << "] " << time_string;
                }
                empty_length_p = empty_length;

                std::string temp_str = strs.str();
                char const* char_type = temp_str.c_str();
                REprintf("\r");
                REprintf("%s", char_type);
            }
            if (_ticks_displayed >= _max_ticks) _finalize_display();
        }
    }

    void flush_console() {
        #if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
            R_FlushConsole();
        #endif
    }

protected:
    int  _max_ticks;
    int  empty_length_p;
    bool _finalized;
    bool _timer_flag;
    time_t start, end;
    int  _ticks_displayed;

    void _finalize_display() {
        if (_finalized) return;
        REprintf("\n");
        _finalized = true;
    }

    std::string _time_to_string(double seconds, float progress) {
        int time = (int) seconds;
        int hour = 0, min = 0, sec = 0;
        hour = time / 3600;
        time = time % 3600;
        min  = time / 60;
        time = time % 60;
        sec  = time;
        std::stringstream time_strs;
        time_strs << (progress < 1.0 ? "TimeLeft: " : "RunTime: ");
        if (hour != 0)              time_strs << hour << "h";
        if (hour != 0 || min != 0)  time_strs << min  << "m";
        time_strs << sec << "s";
        return time_strs.str();
    }

    int _compute_nb_ticks(float progress) {
        return int(progress * _max_ticks);
    }

    virtual std::string _construct_ticks_display_string(int nb) = 0;
};

// ---------------------------------------------------------------------------
// MinimalProgressBar_plus — arrow-style progress bar  [>>>----->      ] ...
// ---------------------------------------------------------------------------
class MinimalProgressBar_plus : public MinimalProgressBarBase {
public:
    MinimalProgressBar_plus() : MinimalProgressBarBase(45) {}
    ~MinimalProgressBar_plus() {}

    void display() override { flush_console(); }

protected:
    std::string _construct_ticks_display_string(int nb) override {
        std::stringstream ticks_strs;
        for (int i = 1; i <= _max_ticks; ++i) {
            if (i < 4) {
                ticks_strs << ">";
            } else if (i < nb) {
                ticks_strs << "-";
            } else if (i == nb) {
                ticks_strs << ">";
            } else {
                ticks_strs << " ";
            }
        }
        return ticks_strs.str();
    }
};

// ---------------------------------------------------------------------------
// MinimalProgressBar_perc — asterisk-style bar with percentage ruler header
// ---------------------------------------------------------------------------
class MinimalProgressBar_perc : public MinimalProgressBarBase {
public:
    MinimalProgressBar_perc() : MinimalProgressBarBase(49) {}
    ~MinimalProgressBar_perc() {}

    void display() override {
        REprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
        REprintf("[----|----|----|----|----|----|----|----|----|----|\n");
        flush_console();
    }

protected:
    std::string _construct_ticks_display_string(int nb) override {
        std::stringstream ticks_strs;
        for (int i = 1; i <= _max_ticks; ++i) {
            if (i <= nb) {
                ticks_strs << "*";
            } else {
                ticks_strs << " ";
            }
        }
        return ticks_strs.str();
    }
};

// ---------------------------------------------------------------------------
// MinimalProgressBar — simple inline percentage display (no ticks)
// ---------------------------------------------------------------------------
class MinimalProgressBar : public ProgressBar {
public:
    MinimalProgressBar()  { _finalized = false; }
    ~MinimalProgressBar() {}

    void display() {}

    void update(float progress) {
        if (_finalized) return;
        REprintf("\r");
        REprintf("Calculating in process...(finished %.2f%%)", progress * 100);
    }

    void end_display() {
        if (_finalized) return;
        REprintf("\r");
        REprintf("Calculating in process...(finished 100.00%%)");
        REprintf("\n");
        _finalized = true;
    }

private:
    bool _finalized;
};

// ---------------------------------------------------------------------------
// fill_geno_buffer — load a batch of cnt markers from a BigMatrix into buf.
//
// buf layout controlled by ind_in_rows:
//   true  → buf is n×cnt, buf(k,l) = genotype of individual k at marker (i_marker+l)
//   false → buf is cnt×n, buf(l,k) = genotype of individual k at marker (i_marker+l)
//
// Empty geno_ind   → use direct row/col index k.
// Empty marker_ind → use sequential index (i_marker + l).
//
// NOTE: Does NOT handle NA values. Callers requiring NA-aware buffer filling
// (e.g. BigRowMean) must use specialized logic.
// ---------------------------------------------------------------------------
template <typename T>
void fill_geno_buffer(arma::mat&         buf,
                      MatrixAccessor<T>& mat,
                      int                cnt,
                      int                i_marker,
                      const arma::uvec&  geno_ind,
                      const arma::uvec&  marker_ind,
                      bool               marker_bycol,
                      bool               ind_in_rows) {
    int  n            = ind_in_rows ? (int)buf.n_rows : (int)buf.n_cols;
    bool has_geno_ind = !geno_ind.is_empty();
    bool has_mrk_ind  = !marker_ind.is_empty();

    if (marker_bycol) {
        if (ind_in_rows) {
            #pragma omp parallel for
            for (int l = 0; l < cnt; l++) {
                int col = has_mrk_ind ? (int)marker_ind[i_marker + l] : (i_marker + l);
                for (int k = 0; k < n; k++)
                    buf(k, l) = (double)mat[col][has_geno_ind ? (int)geno_ind[k] : k];
            }
        } else {
            #pragma omp parallel for
            for (int l = 0; l < cnt; l++) {
                int col = has_mrk_ind ? (int)marker_ind[i_marker + l] : (i_marker + l);
                for (int k = 0; k < n; k++)
                    buf(l, k) = (double)mat[col][has_geno_ind ? (int)geno_ind[k] : k];
            }
        }
    } else {
        if (ind_in_rows) {
            #pragma omp parallel for
            for (int k = 0; k < n; k++) {
                int row = has_geno_ind ? (int)geno_ind[k] : k;
                for (int l = 0; l < cnt; l++)
                    buf(k, l) = (double)mat[row][has_mrk_ind ? (int)marker_ind[i_marker + l] : (i_marker + l)];
            }
        } else {
            #pragma omp parallel for
            for (int k = 0; k < n; k++) {
                int row = has_geno_ind ? (int)geno_ind[k] : k;
                for (int l = 0; l < cnt; l++)
                    buf(l, k) = (double)mat[row][has_mrk_ind ? (int)marker_ind[i_marker + l] : (i_marker + l)];
            }
        }
    }
}

// Dispatch a bigmatrix function call based on element type.
// Requires xpMat to be the first argument to FUNC<T>; trailing args forwarded via __VA_ARGS__.
// Usage: DISPATCH_MATRIX_TYPE(myFunc, xpMat, arg1, arg2, ...)
// Note: __VA_ARGS__ is intentionally never empty at any call site, so no ## GNU extension needed.
#define DISPATCH_MATRIX_TYPE(FUNC, xpMat, ...)                              \
    switch((xpMat)->matrix_type()) {                                         \
    case 1: return FUNC<char>  (xpMat, __VA_ARGS__);                       \
    case 2: return FUNC<short> (xpMat, __VA_ARGS__);                       \
    case 4: return FUNC<int>   (xpMat, __VA_ARGS__);                       \
    case 8: return FUNC<double>(xpMat, __VA_ARGS__);                       \
    default: throw Rcpp::exception("unsupported bigmatrix type");           \
    }

#endif  // RMVP_H_
