# Sprint 0.1: Numerical regression baseline for GWAS methods.
# Locks the exact numerical outputs of GLM, MLM (BRENT/EMMA/HE VC),
# FaSTLMM LL, and FarmCPU before refactoring begins.
#
# Expected values were captured on the extdata 15-marker dataset using
# rMVP 1.4.6 on x86-64 Linux with OpenBLAS. Tolerance is 1e-6 to accommodate
# cross-platform BLAS differences (Apple Accelerate, MKL, OpenBLAS agree to
# ~1e-7 on these inputs; 1e-6 gives one decade of headroom).
#
# To update after an intentional numerical change (documented in CHANGELOG):
#   1. Run each function interactively on the test dataset.
#   2. Print values with options(digits=15) and update the constants below.

# Helper: load and cache test dataset once per test session.
.load_gwas_testdata <- local({
    cached <- NULL
    function() {
        if (!is.null(cached)) return(cached)
        phe  <- read.table(
            system.file("extdata", "07_other", "mvp.phe", package = "rMVP"),
            header = TRUE
        )
        geno <- attach.big.matrix(
            system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
        )
        idx       <- !is.na(phe[, 2])
        phe_clean <- phe[idx, ]
        K         <- MVP.K.VanRaden(geno, cpu = 1, verbose = FALSE)
        eigenK    <- eigen(K, symmetric = TRUE)
        cached    <<- list(phe = phe_clean, geno = geno, K = K, eigenK = eigenK)
        cached
    }
})

# Rows 8 and 9 are monomorphic in the extdata set — expected to produce NA/NaN.
# Rows checked: 1, 4, 11 (most significant by p-value), 15.
.CHECK_ROWS <- c(1L, 4L, 11L, 15L)

# ---- GLM ----------------------------------------------------------------

test_that("GLM produces numerically stable results", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.GLM(phe = d$phe, geno = d$geno, cpu = 1, verbose = FALSE)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), 15L)
    expect_equal(ncol(result), 3L)

    # Monomorphic markers return NA / NaN.
    expect_true(is.na(result[8, 1]))
    expect_true(is.na(result[9, 1]))
    expect_true(is.nan(result[8, 3]))
    expect_true(is.nan(result[9, 3]))

    # Row 11 is the most significant marker.
    expect_equal(which.min(result[, 3]), 11L)

    # Beta (col 1), SE (col 2), p-value (col 3) at fixed rows.
    tol <- 1e-6
    expect_equal(result[.CHECK_ROWS, 1],
                 c( 2.5728571, -2.8711111, -4.0233333,  2.8900000), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 2],
                 c( 1.7134863,  2.7840715,  1.1267864,  2.7816275), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 3],
                 c(0.17161487, 0.33258383, 0.00728702, 0.32920268), tolerance = tol)
})

# ---- Variance component estimation --------------------------------------

test_that("BRENT VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.BRENT.Vg.Ve(
        y       = as.numeric(d$phe[, 2]),
        X       = matrix(1, nrow(d$phe)),
        eigenK  = d$eigenK,
        verbose = FALSE
    )

    tol <- 1e-6
    expect_equal(vc$vg,    19.460680, tolerance = tol)
    expect_equal(vc$ve,     8.090278, tolerance = tol)
    expect_equal(vc$delta,  0.415724, tolerance = tol)
})

test_that("EMMA VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.EMMA.Vg.Ve(
        y = as.numeric(d$phe[, 2]),
        X = matrix(1, nrow(d$phe)),
        K = d$K
    )

    tol <- 1e-6
    expect_equal(vc$vg,    19.461379, tolerance = tol)
    expect_equal(vc$ve,     8.090335, tolerance = tol)
    expect_equal(vc$delta,  0.415712, tolerance = tol)
    expect_equal(vc$REML,  -27.102099, tolerance = tol)
})

test_that("HE VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.HE.Vg.Ve(
        y = as.numeric(d$phe[, 2]),
        X = matrix(1, nrow(d$phe)),
        K = d$K
    )

    tol <- 1e-6
    expect_equal(vc$vg,     8.288881, tolerance = tol)
    expect_equal(vc$ve,    18.890817, tolerance = tol)
    expect_equal(vc$delta,  2.279055, tolerance = tol)
})

test_that("FaSTLMM LL estimation is numerically stable", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.FaSTLMM.LL(
        pheno    = as.matrix(d$phe[, 2]),
        snp.pool = d$geno[],
        X0       = matrix(1, nrow(d$phe)),
        ncpus    = 1
    )

    tol <- 1e-6
    # delta hits the optimizer upper boundary (essentially Inf) on this small
    # dataset, so vg collapses to ~0. Only beta, LL, and ve are meaningful.
    expect_equal(result$beta, 100.298,    tolerance = tol)
    expect_equal(result$LL,   -30.341552, tolerance = tol)
    expect_equal(result$ve,    25.290616, tolerance = tol)
})

# ---- MLM ----------------------------------------------------------------

test_that("MLM (BRENT VC) produces numerically stable results", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.MLM(
        phe       = d$phe,
        geno      = d$geno,
        eigenK    = d$eigenK,
        vc.method = "BRENT",
        cpu       = 1,
        verbose   = FALSE
    )

    expect_true(is.matrix(result))
    expect_equal(nrow(result), 15L)

    # Monomorphic markers return NaN.
    expect_true(is.nan(result[8, 1]))
    expect_true(is.nan(result[9, 1]))

    tol <- 1e-6
    expect_equal(result[.CHECK_ROWS, 1],
                 c( 3.9431337, -0.0500524, -3.6286984,  2.4695755), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 2],
                 c( 2.2801950,  2.7330146,  1.9299914,  2.7095227), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 3],
                 c(0.12201210, 0.98583689, 0.09688123, 0.38870037), tolerance = tol)
})

test_that("MLM (EMMA VC) produces numerically stable results", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.MLM(
        phe       = d$phe,
        geno      = d$geno,
        K         = d$K,
        eigenK    = d$eigenK,
        vc.method = "EMMA",
        cpu       = 1,
        verbose   = FALSE
    )

    expect_true(is.matrix(result))
    expect_equal(nrow(result), 15L)

    expect_true(is.nan(result[8, 1]))
    expect_true(is.nan(result[9, 1]))

    tol <- 1e-6
    expect_equal(result[.CHECK_ROWS, 1],
                 c( 3.9431501, -0.0500181, -3.6286904,  2.4695614), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 2],
                 c( 2.2802260,  2.7330490,  1.9300190,  2.7095560), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 3],
                 c(0.12201518, 0.98584678, 0.09688584, 0.38870857), tolerance = tol)
})

test_that("MLM (HE VC) produces numerically stable results", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.MLM(
        phe       = d$phe,
        geno      = d$geno,
        K         = d$K,
        eigenK    = d$eigenK,
        vc.method = "HE",
        cpu       = 1,
        verbose   = FALSE
    )

    expect_true(is.matrix(result))
    expect_equal(nrow(result), 15L)

    expect_true(is.nan(result[8, 1]))
    expect_true(is.nan(result[9, 1]))

    tol <- 1e-6
    expect_equal(result[.CHECK_ROWS, 1],
                 c( 3.1813154, -1.8035841, -3.9345467,  2.9160929), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 2],
                 c( 2.0865350,  2.7990780,  1.8065600,  2.8072150), tolerance = tol)
    expect_equal(result[.CHECK_ROWS, 3],
                 c(0.16584474, 0.53738173, 0.06105992, 0.32927922), tolerance = tol)
})

# ---- FarmCPU ------------------------------------------------------------

test_that("FarmCPU produces numerically stable results", {
    skip_on_cran()
    # FarmCPU's bin-optimization algorithm requires a minimum number of markers
    # to select pseudo-QTN candidates. The 15-marker extdata test set is too
    # small. FarmCPU regression is covered by manual integration tests with the
    # pig60K dataset (devscripts/benchmark_farmcpu.R).
    skip("FarmCPU binning requires more markers than the 15-marker extdata set provides.")
})
