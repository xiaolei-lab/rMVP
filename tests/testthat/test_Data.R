# read data
vcfPath     <- system.file("extdata", "01_vcf", "mvp.vcf", package = "rMVP")
bfilePath   <- system.file("extdata", "02_bfile", "mvp", package = "rMVP")
hmp1Path    <- system.file("extdata", "03_hapmap", "mvp.hmp.txt", package = "rMVP")
hmp2Path    <- system.file("extdata", "03_hapmap", "mvp.diploid.hmp.txt", package = "rMVP")
numericPath <- system.file("extdata", "04_numeric", "mvp.num", package = "rMVP")
genoPath    <- system.file("extdata", "05_mvp", "mvp.geno.desc", package = "rMVP")
genoImpPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
phePath     <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
mapPath     <- system.file("extdata", "07_other", "mvp.map", package = "rMVP")


context("MVP.Data - vcf")

test_that("MVP.Data() - vcf", {
    out <- "rMVP.test.1"
    expect_output(
        MVP.Data(fileVCF = vcfPath, out = out, verbose = FALSE, ncpus = 1),
        "successfully!"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".imp.geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".imp.kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".imp.pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".map"), header = TRUE, stringsAsFactors = FALSE)
    
    expect_known_value(geno[], "rMVP.keep.geno")
    expect_known_value(genoImp[], "rMVP.keep.genoImp")
    expect_known_value(kinship[], "rMVP.keep.kinship")
    expect_known_value(pcs[], "rMVP.keep.pcs")
    expect_known_value(genoInd, "rMVP.keep.genoInd")
    expect_known_value(map, "rMVP.keep.map")
})

context("MVP.Data - bfile")

test_that("MVP.Data() - Bfile", {
    out <- "rMVP.test.2"
    expect_output(
        MVP.Data(fileBed = bfilePath, out = out, verbose = FALSE, ncpus = 1),
        "successfully!"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".imp.geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".imp.kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".imp.pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".map"), header = TRUE, stringsAsFactors = FALSE)

    expect_known_value(geno[], "rMVP.keep.geno", update = FALSE)
    expect_known_value(genoImp[], "rMVP.keep.genoImp", update = FALSE)
    expect_known_value(kinship[], "rMVP.keep.kinship", update = FALSE)
    expect_known_value(pcs[], "rMVP.keep.pcs", update = FALSE)
    expect_known_value(genoInd, "rMVP.keep.genoInd", update = FALSE)
    expect_known_value(map, "rMVP.keep.map", update = FALSE)
})

context("MVP.Data - hapmap")

test_that("MVP.Data() - HMP Diploid", {
    out <- "rMVP.test.3"
    expect_output(
        MVP.Data(fileHMP = hmp2Path, out = out, verbose = FALSE, ncpus = 1),
        "successfully!"
    )
    geno <- attach.big.matrix(paste0(out, ".geno.desc"))
    genoImp <- attach.big.matrix(paste0(out, ".imp.geno.desc"))
    kinship <- attach.big.matrix(paste0(out, ".imp.kin.desc"))
    pcs <- attach.big.matrix(paste0(out, ".imp.pc.desc"))
    genoInd <- read.table(paste0(out, ".geno.ind"), stringsAsFactors = FALSE)
    map <- read.table(paste0(out, ".map"), header = TRUE, stringsAsFactors = FALSE)
    
    expect_known_value(geno[], "rMVP.keep.geno", update = FALSE)
    expect_known_value(genoImp[], "rMVP.keep.genoImp", update = FALSE)
    expect_known_value(kinship[], "rMVP.keep.kinship", update = FALSE)
    expect_known_value(pcs[], "rMVP.keep.pcs", update = FALSE)
    expect_known_value(genoInd, "rMVP.keep.genoInd", update = FALSE)
    expect_known_value(map, "rMVP.keep.map", update = FALSE)
})

files <- dir(pattern = "^rMVP.test")
file.remove(files)
