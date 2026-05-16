# Package-level constants for rMVP statistical methods.
# Naming convention: .<METHOD>_<DESCRIPTION> (dot-prefix keeps them non-exported)

# FaSTLMM log-likelihood optimization
.FASTLMM_SVD_THRESHOLD        <- 1e-8   # singular value cutoff in SVD
.FASTLMM_DELTA_EXP_START      <- -5     # log-delta grid lower bound
.FASTLMM_DELTA_EXP_END        <- 5      # log-delta grid upper bound
.FASTLMM_DELTA_EXP_STEP       <- 0.1    # log-delta grid step size
.FASTLMM_DELTA_EXP_DEGENERATE <- 100    # sentinel: collapse grid to single point when SNP pool has a constant column
.FASTLMM_2PI                  <- 2 * pi # NOTE: original code had 2 * 3.14 (bug); corrected to 2 * pi

# EMMA variance component estimation (defaults match original function signature)
.EMMA_NGRIDS <- 100
.EMMA_LLIM   <- -10
.EMMA_ULIM   <- 10
.EMMA_ESP    <- 1e-10

# BRENT variance component estimation
.BRENT_MIN_EIGENVALUE <- 1e-6   # eigenvalue floor to avoid division by near-zero

# HE regression
.HE_DEFAULT_LOG_SIGMA2 <- log(0.1)   # log-sigma2 fallback when CalcVChe returns non-positive value

# FarmCPU
.FARMCPU_LD_THRESHOLD <- 0.7    # LD threshold for pseudo-QTN deduplication
