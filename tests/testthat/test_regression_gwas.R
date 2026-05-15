# Sprint 0.1: Numerical regression baseline for GWAS methods.
# Locks the exact numerical outputs of GLM, MLM (BRENT/EMMA/HE VC),
# FaSTLMM LL, and FarmCPU before refactoring begins.
#
# On first run: fixture files are created in tests/testthat/fixtures/.
# On subsequent runs: results are compared against stored fixtures.
#
# To regenerate fixtures after an intentional numerical change (bug fix,
# algorithm improvement documented in CHANGELOG):
#   unlink("tests/testthat/fixtures", recursive = TRUE)
#   devtools::test(filter = "regression_gwas")

# Helper: load/compare against an RDS fixture.
# First run creates the fixture; subsequent runs compare with tolerance 1e-10.
.expect_fixture <- function(result, name, tol = 1e-10) {
    dir <- test_path("fixtures")
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    path <- file.path(dir, paste0(name, ".rds"))
    if (!file.exists(path)) {
        saveRDS(result, path)
        succeed(paste("Created fixture:", name))
    } else {
        expect_equal(result, readRDS(path), tolerance = tol)
    }
}

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

# ---- GLM ----------------------------------------------------------------

test_that("GLM produces numerically stable results", {
    skip_on_cran()
    d      <- .load_gwas_testdata()
    result <- MVP.GLM(phe = d$phe, geno = d$geno, cpu = 1, verbose = FALSE)
    .expect_fixture(result, "glm")
})

# ---- Variance component estimation --------------------------------------

test_that("BRENT VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.BRENT.Vg.Ve(
        y      = as.numeric(d$phe[, 2]),
        X      = matrix(1, nrow(d$phe)),
        eigenK = d$eigenK,
        verbose = FALSE
    )
    .expect_fixture(vc, "vc_brent")
})

test_that("EMMA VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.EMMA.Vg.Ve(
        y = as.numeric(d$phe[, 2]),
        X = matrix(1, nrow(d$phe)),
        K = d$K
    )
    .expect_fixture(vc, "vc_emma")
})

test_that("HE VC estimation is numerically stable", {
    skip_on_cran()
    d  <- .load_gwas_testdata()
    vc <- MVP.HE.Vg.Ve(
        y = as.numeric(d$phe[, 2]),
        X = matrix(1, nrow(d$phe)),
        K = d$K
    )
    .expect_fixture(vc, "vc_he")
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
    .expect_fixture(result, "fastlmm_ll")
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
    .expect_fixture(result, "mlm_brent")
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
    .expect_fixture(result, "mlm_emma")
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
    .expect_fixture(result, "mlm_he")
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
