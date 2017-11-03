MVP.PCA <-
function(M, perc=1, pcs.keep=5, memo=NULL){
##########################################################################################################
# Object: Principal Component Analysis
# Input:
# M: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
# perc: percentage of total SNPs selected for PCA
# pcs.keep: maximum number of PCs
# memo: a marker on output file name
#
# Output:
# PCs
# Authors: Xiaolei Liu and Lilin Yin
# Build date: Dec 14, 2016
# Last update: Dec 14, 2016
##########################################################################################################
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    linux <- R.ver == 'Linux'
    mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")

    if(r.open && mac){
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    #Data Check
    if(is.null(M)){
        stop("There is no genotype data!")
    }

    if(pcs.keep < 0){
        pcs.keep = 0
    }

    if(pcs.keep == 0){
        return(NULL)
    }

    if(perc < 0 | perc > 1){
        print("perc value should less than 1 and bigger than 0, MVP reset it to 1")
        perc <- 1
    }
    print("Principal Component Analysis Start...")
    if(perc < 1){
        if(wind){
            set.seed(123456)
            m <- nrow(M)
            big.geno <- M[sample(1:m, round(m * perc)), ]
        }else{
            big.geno <- thin(bigMat=M, keep = perc)
        }
    }else{
        if(wind){
            big.geno <- t(as.matrix(M))
        }else{
            big.geno <- M
        }
    }
    #Calculate PCs by using bigPCA package
    if(wind){
        PCs <- prcomp(big.geno)$x[, 1:pcs.keep]
    }else{
        PCs <- big.PCA(bigMat=big.geno,pcs.to.keep=pcs.keep)$PCs
    }
    return(list(PCs=PCs))
}#end of MVP.PCA function
