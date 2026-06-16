# read data
vcfPath     <- system.file("extdata", "01_vcf", "mvp.vcf", package = "rMVP")
bfilePath   <- file.path(system.file("extdata", "02_bfile", package = "rMVP"), "mvp")
hmp1Path    <- system.file("extdata", "03_hapmap", "mvp.hmp.txt", package = "rMVP")
numericPath <- system.file("extdata", "04_numeric", "mvp.num", package = "rMVP")
genoPath    <- system.file("extdata", "05_mvp", "mvp.geno.desc", package = "rMVP")
genoImpPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
phePath     <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
mapPath     <- system.file("extdata", "07_other", "mvp.map", package = "rMVP")


context("MVP.Data - vcf")

test_that("MVP.Data() - vcf", {
    skip_on_cran()
    
    out <- "rMVP.test.1"
    expect_output(
        MVP.Data(fileVCF = vcfPath, out = out, fileKin = TRUE, filePC = TRUE, verbose = FALSE, ncpus = 2, maxLine = 1e3),
        "done"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".geno.map"), header = TRUE, stringsAsFactors = FALSE)
    
    expect_known_value(geno[], "rMVP.keep.geno")
    expect_known_value(genoImp[], "rMVP.keep.genoImp")
    expect_known_value(kinship[], "rMVP.keep.kinship")
    expect_known_value(pcs[], "rMVP.keep.pcs")
    expect_known_value(genoInd, "rMVP.keep.genoInd")
    expect_known_value(map, "rMVP.keep.map")
})

context("MVP.Data - bfile")

test_that("MVP.Data() - Bfile", {
    skip_on_cran()
    
    out <- "rMVP.test.2"
    expect_output(
        MVP.Data(fileBed = bfilePath, out = out, fileKin = TRUE, filePC = TRUE, verbose = FALSE, ncpus = 2, maxLine = 1e3),
        "done"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".geno.map"), header = TRUE, stringsAsFactors = FALSE)

    expect_known_value(geno[], "rMVP.keep.geno", update = FALSE)
    expect_known_value(genoImp[], "rMVP.keep.genoImp", update = FALSE)
    expect_known_value(kinship[], "rMVP.keep.kinship", update = FALSE)
    expect_known_value(pcs[], "rMVP.keep.pcs", update = FALSE)
    expect_known_value(genoInd, "rMVP.keep.genoInd", update = FALSE)
    expect_known_value(map, "rMVP.keep.map", update = FALSE)
})

context("MVP.Data - hapmap")

test_that("MVP.Data() - HMP", {
    skip_on_cran()
    
    out <- "rMVP.test.3"
    expect_output(
        MVP.Data(fileHMP = hmp1Path, out = out, fileKin = TRUE, filePC = TRUE, verbose = FALSE, ncpus = 2, maxLine = 1e3),
        "done"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".geno.map"), header = TRUE, stringsAsFactors = FALSE)
    
    expect_known_value(geno[], "rMVP.keep.geno", update = FALSE)
    expect_known_value(genoImp[], "rMVP.keep.genoImp", update = FALSE)
    expect_known_value(kinship[], "rMVP.keep.kinship", update = FALSE)
    expect_known_value(pcs[], "rMVP.keep.pcs", update = FALSE)
    expect_known_value(genoInd, "rMVP.keep.genoInd", update = FALSE)
    expect_known_value(map, "rMVP.keep.map", update = FALSE)
})

context("MVP.Data - numeric")

test_that("MVP.Data.Numeric2MVP() - Marker-by-Individual format (Bug #116 fix)", {
    skip_on_cran()
    
    # Test for the SetRows.bm error fix in MVP.Data.Numeric2MVP
    # The example numeric file has 15 markers (rows) x 10 individuals (columns)
    # This tests the Marker-by-Individual format detection and assignment logic
    
    out <- "rMVP.test.numeric"
    mapPath <- system.file("extdata", "04_numeric", "mvp.map", package = "rMVP")
    numericPath <- system.file("extdata", "04_numeric", "mvp.num", package = "rMVP")
    
    # This should NOT raise "Illegal row index usage in extraction" error
    expect_silent(
        MVP.Data.Numeric2MVP(numericPath, mapPath, out = out, verbose = FALSE)
    )
    
    # Verify the output files were created
    expect_true(file.exists(paste0(out, ".geno.desc")))
    expect_true(file.exists(paste0(out, ".geno.bin")))
    expect_true(file.exists(paste0(out, ".geno.map")))
    expect_true(file.exists(paste0(out, ".geno.ind")))
    
    # Verify dimensions: 10 individuals x 15 markers
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    expect_equal(nrow(geno), 10L, label = "number of individuals")
    expect_equal(ncol(geno), 15L, label = "number of markers")
    
    # Verify data integrity - check first few values match input
    # The first individual should have genotypes: 0, 2, 2, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0
    expect_equal(as.numeric(geno[1, ]), 
                 c(0, 2, 2, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0))
})

files <- dir(pattern = "^rMVP.test")
file.remove(files)
files <- dir(pattern = "*.log")
file.remove(files)