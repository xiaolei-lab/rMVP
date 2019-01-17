# Data pre-processing module
# 
# Copyright (C) 2016-2018 by Xiaolei Lab
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Principal Component Analysis
#' 
#' Build date: Dec 14, 2016
#' Last update: Oct 29, 2018
#' @author Xiaolei Liu, Lilin Yin and Haohao Zhang
#' 
#' @param M Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param perc percentage of total SNPs selected for calculate principal components
#' @param pcs.keep maximum number of PCs for output
#' @param memo a marker on output file name
#'
#' @return
#' Output: PCs - a n * npc matrix of top number of PCs, n is population size and npc is @param pcs.keep 
#' 
#' @export
#'
#' @examples
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' pca <- MVP.PCA(M=genotype)
#' str(pca)
MVP.PCA <-
function(M, perc=1, pcs.keep=5, memo=NULL){
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    linux <- R.ver == 'Linux'
    mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent = TRUE),"try-error")

    if (r.open && mac) {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    #Data Check
    if (is.null(M)) {
        stop("There is no genotype data!")
    }

    if (pcs.keep < 0) {
        pcs.keep = 0
    }

    if (pcs.keep == 0) {
        return(NULL)
    }

    if (perc < 0 | perc > 1) {
        print("perc value should less than 1 and bigger than 0, MVP reset it to 1")
        perc <- 1
    }
    print("Principal Component Analysis Start...")
    if (perc < 1) {
        # set.seed(123456)
        m <- nrow(M)
        big.geno <- t(M[sample(1:m, round(m * perc)), ])
    }else{
        big.geno <- t(as.matrix(M))
    }
    #Calculate PCs by using bigPCA package

    PCs <- prcomp(big.geno)$x[, 1:pcs.keep]

    return(list(PCs=PCs))
}#end of MVP.PCA function


# Copy From https://github.com/nicholasjcooper/bigpca/blob/master/bigpca.R
#' PCA/Singular Value Decomposition for big.matrix
#' 
#' At the time of submission there was no existing native method to conduct principle components
#' analysis (PCA) on big.matrix objects. This function allows singular value decomposition (SVD) of
#' very large matrices, very efficiently, versus the default method. The major speed advantages
#' occur when the 'bigalgebra' package is installed, and when the argument for this function
#' 'SVD'=TRUE. Regular PCA can be conducted using SVD=FALSE but it will be slower and the maximum
#' matrix size able to produce a result, given memory limits, will be smaller. SVD is not exactly
#' the same as PCA, but from my testing the components produced will correlate R>.9 with components
#' of PCA on the same matrix. This function is not completely native to big.matrix objects so 
#' there is one step where the matrix submitted needs to be loaded into memory, so if your 
#' big.matrix object is larger than the allowed size of a standard R-matrix [[which is roughly 3GB; 
#' you can check using NCmisc::estimate.memory()], then this function will fail unless you set the 
#' option 'thin' to a percentage that, multiplied by the original matrix memory-load, is under 3GB.
#' For large matrices in my applications, components produced with thinning are still highly 
#' correlated with components produced using the full dataset. For a breakdown of thinning methods,
#' see the description for the function thin() for more  information. Even with medium sized
#' matrices, for instance 15,000 x 50,000 in size, this function is orders of magnitude faster
#' than the standard R PCA functions, usually running in a matter of minutes, rather than hours
#' or days in examples that I have tested, due to much better handling of memory for internal
#' transpose and eigen operations by using the 'bigmemory' architecture.
#'
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param dir directory containing the filebacked.big.matrix, and also where the output
#'  file will be saved by default if the save.pcs option is set TRUE. 
#' @param pcs.to.keep integer, number of principle components to keep. Singular Value Decomposition
#'  methods are able to run faster if they don't need to calculate every single PC for a large
#'  matrix. Default is to calculate only the first 50; in practice even fewer than this are generally
#'  used directly. Apart from reducing processing time, this can also reduce storage/RAM burden for 
#'  the resulting matrix. Set to NA, or a number >= min(dim(bigMat)) in order to keep all PCs.
#' @param thin decimal, percentage of the original number of rows you want to thin the matrix to.
#'  see function thin() for details of how this can be done, pass arguments to thin() using ...
#'  Even though this PCA function uses mainly 'big.matrix' native methods, there is a step where the
#'  matrix must be stored fully in memory, so this limits the size of what matrix can be processed,
#'  depending on RAM limits. If you want to conduct PCA/SVD on a matrix larger than RAM you can thin
#'  the matrix to a percentage of the original size. Usually such a large matrix will contain correlated
#'  measures and so the exclusion of some data-rows (variables) will have only a small impact on the
#'  resulting principle components. In some applications tested using this function, using only 5% 
#'  of 200,000 variables a PCA gave extremely similar results to using the full dataset.
#' @param SVD logical, whether to use a Singular Value Decomposition method or a PCA method. The 
#'  eigenvalues and eigenvectors of each alternative will be highly correlated so for most applications,
#'  such as PC-correction, this shouldn't make much difference to the result. However, using SVD=TRUE
#'  can provide significant advantages in speed, or even allow a solution on a matrix that would be
#'  to large to practically compute full PCA. Note that only in SVD mode, and with the bigalgebra
#'  package installed will the full speed advantage of this function be utilised.
#' @param LAP logical, whether to use La.svd() instead of svd() when SVD=TRUE, see base:svd for more info.
#' @param center whether to 'centre' the matrix rows by subtracting rowMeans() before conducting the PCA. This is usually
#'  advisable, although you may wish to skip this if the matrix is already centred to save extra processing.
#'  unlike prcomp there is no option to standardize or use the correlation matrix, if this is desired, please
#'  standardize the bigMat object before running this function. Alternatively, 'center' can be a vector of length
#'  nrow(bigMat) which gives an offset to correct each row by.
#' @param save.pcs whether to save the principle component matrix and eigenvalues to a binary file with name pcs.fn
#' @param pcs.fn name of the binary when save.pcs=TRUE
#' @param return.loadings logical, whether to return the 'loadings' (or the other singular vector when SVD=TRUE); could result in a speed decrease
#' @param verbose whether to display detailed progress of the PCA
#' @param use.bigalgebra logical, whether to use the bigalgebra package for algebraic operations. For large
#'  datasets bigalgebra should provide a substantial speedup, and also facilitates use of larger matrices. This 
#'  relies on having bigalgebra installed and loaded, which requires some manual configuration as bigalgebra
#'  is not currently available on CRAN, but only SVN and RForge. See svn.bigalgebra.install() or big.algebra.install.help()
#'  Default is to use bigalgebra if available (TRUE), but setting this FALSE prevents the check for bigalgebra which would be
#'  cleaner if you know that you don't have it installed.
#' @param ... if thin is TRUE, then these should be any additional arguments for thin(), e.g, 'pref', 'keep', 'how', etc.i
#' @param delete.existing logical, whether to automatically delete filebacked matrices (if they exist) 
#' before rewriting. This is because of an update since 20th October 2015 where bigmemory won't allow
#' overwrite of an existing filebacked matrix. If you wish to set this always TRUE or FALSE, use
#'  options(deleteFileBacked)
#' @return A list of principle components/singular vectors (may be incomplete depending on options selected), and of
#'  the eigenvalues/singular values.
#' @export
#' @seealso \code{\link{get.big.matrix}}, \code{\link{PC.correct}}
#' @author Nicholas Cooper
#' @examples 
## PRELIMINARY EXAMPLES: demonstration of PCA versus SVD ##
#' # create an example matrix and its transpose
#' min.dim <- 200; nvar <- 500; subset.size <- 50
#' mat <- matrix(rnorm(min.dim*nvar),ncol=min.dim) 
#' prv.large(mat)
#' t.mat <- t(mat)
#' # create two alternative covariance matrices
#' MMs <- t.mat %*% mat
#' MsM <- mat %*% t.mat
#' # run singular value decomposition
#' pca <- svd(mat)   
#' D <- pca$d # singular values (=sqrt(eigenvalues))
#' V <- pca$v # right singular vector
#' U <- pca$u # left singular vector
#' sig <- mat-mat; diag(sig) <- D; 
#' MMs2 <- V %*% (t(sig) %*% sig) %*% t(V)
#' sig <- t.mat-t.mat; diag(sig) <- D; 
#' MsM2 <- U %*% (sig %*% t(sig)) %*% t(U)
#' # show that the covariance matrices are equal to the functions of 
#' # the left and right singular vectors
#' prv(MMs,MsM); prv(MMs2,MsM2)
#' pr <- princomp(mat) # PCA using eigendecomposition of cov matrix
#' L <- matrix(rep(0,40000),ncol=200); diag(L) <- pr[[1]]^2 # eigenvalues as diag
#' mat2 <- (pr[[2]]) %*% L %*%  solve(pr[[2]]) # = eigenvectors * eigenvalues * inv(eigenvectors)
#' prv.large(cov(mat)); prv.large(mat2) #  == COVmat (may be slight tolerance differences)
#' ## Now demonstrate the correlation between SVD and PCA ##
#' # the right singular vector is highly correlated with the pca loadings:
#' median(abs(diag(cor(V,pr[["loadings"]]))))
#' # the left singular vector is highly correlated with the pca scores (eigenvectors):
#' median(abs(diag(cor(U,pr[["scores"]]))))
#' cor(pr$sdev,D) # the singular values are equivalent to the eigenvalues
#' 
#' ## MAIN EXAMPLES ##
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' if(file.exists("testMyBig.bck")) { unlink(c("testMyBig.bck","testMyBig.dsc")) }
#' bmat <- as.big.matrix(mat,backingfile="testMyBig.bck",
#'  descriptorfile="testMyBig.dsc",  backingpath = getwd())
#' result <- big.PCA(bmat) #,verbose=TRUE)
#' headl(result)
#' # plot the eigenvalues with a linear fit line and elbow placed at 13
#' Eigv <- pca.scree.plot(result$Evalues,M=bmat,elbow=6,printvar=FALSE)
#' rm(bmat) 
#' unlink(c("testMyBig.bck","testMyBig.dsc"))
#' ##  generate some data with reasonable intercorrelations ##
#' mat2 <- sim.cor(500,200,genr=function(n){ (runif(n)/2+.5) })
#' bmat2 <- as.big.matrix(mat2,backingfile="testMyBig2.bck",
#'  descriptorfile="testMyBig2.dsc",  backingpath = getwd())
#' # calculate PCA on decreasing subset size 
#' result2 <- big.PCA(bmat2,thin=FALSE)
#' normal <- result2$PCs; rm(result2)
#' result3 <- big.PCA(bmat2,thin=TRUE,keep=.5, pref="t1")
#' thinned <- result3$PCs; rm(result3)
#' result4 <- big.PCA(bmat2,thin=TRUE,keep=.5, pref="t2", how="cor")
#' corred <- result4$PCs; rm(result4)
#' result5 <- big.PCA(bmat2,thin=TRUE,keep=.5, pref="t3", how="pca")
#' pced <- result5$PCs; rm(result5)
#' result6 <- big.PCA(bmat2,thin=TRUE,keep=.2, pref="t4")
#' thinner <- result6$PCs; rm(result6)
#' ## correlate the resulting PCs with the un-thinned PCs
#' cors.thin.with.orig <- apply(cor(normal,thinned),1,max)
#' cors.corred.with.orig <- apply(cor(normal,corred),1,max)
#' cors.pced.with.orig <- apply(cor(normal,pced),1,max)
#' cors.thinner.with.orig <-apply(cor(normal,thinner),1,max)
#' plot(cors.thin.with.orig,type="l",col="red",ylim=c(0,1))
#' lines(cors.thinner.with.orig,col="orange")
#' lines(cors.corred.with.orig,col="lightblue")
#' lines(cors.pced.with.orig,col="lightgreen")
#' # can see that the first component is highly preserved,
#' # and next components, somewhat preserved; try using different thinning methods
#' rm(bmat2) 
#' unlink(c("testMyBig2.bck","testMyBig2.dsc"))
#' setwd(orig.dir)
# big.PCA <- function(bigMat, dir=getwd(), pcs.to.keep=50, thin=FALSE, SVD=TRUE, LAP=FALSE, center=TRUE,
#                     save.pcs=FALSE, use.bigalgebra=TRUE, pcs.fn="PCsEVsFromPCA.RData", return.loadings=FALSE,
#                     verbose=FALSE, delete.existing=getOption("deleteFileBacked"), ...) {
#     # run principle components analysis on the SNP subset of the LRR snp x sample matrix
#     # various methods to choose from with pro/cons of speed/memory, etc.
#     #  must use SNP-subset to avoid LD, destroying main effects, +avoid huge memory requirements
#     #dir <- validate.dir.for(dir,c("big","pc"))
#     if (exists("validate.dir.for",mode = "function")) {
#         ## plumbCNV specific code ##
#         dir <- do.call("validate.dir.for",list(dir = dir,elements = c("big","pc"), warn = FALSE))  
#     } else {
#         # otherwise
#         dir <- list(big = dir, pc = dir)
#     }
#     #must.use.package(c("irlba"),TRUE)
#     if (thin) {
#         if (verbose) {  prv.big.matrix(bigMat) }
#         bigMat <- thin(bigMat, dir = dir,...)
#         if (verbose) {  prv.big.matrix(bigMat) }
#     } 
#     pcaMat <- get.big.matrix(bigMat,dir)
#     #print(dim(pcaMat))
#     if (verbose & !thin) { prv.big.matrix(pcaMat,name = "pcaMat") }
#     est.mem <- estimate.memory(pcaMat)
#     if (est.mem > 1) {
#         cat(" estimated memory required for", nrow(pcaMat), "x", ncol(pcaMat), "matrix:", round(est.mem,2),
#             "GB. If this exceeds available,\n  then expect PCA to take a long time or fail!\n")
#     }
#     #print(packages.loaded())
#     subMat <- pcaMat[1:nrow(pcaMat),1:ncol(pcaMat)] # must convert bigmatrix to plain matrix here, no pca yet takes a bigmatrix
#     rm(pcaMat)
#     # center using row means
#     if (length(center) > 1) {
#         if (length(center) == nrow(subMat)) {
#             CM <- center
#             center <- TRUE
#             rm.sub <- function(X) { 
#                 mmm <-  matrix(rep(CM, times = ncol(subMat)), ncol = ncol(subMat),byrow = FALSE)
#                 prv(mmm)
#                 return(mmm)
#             }
#         } else {
#             rm.sub <- rowMeans ; warning("center vector didn't match number of rows of 'bigMat', data left uncentered")
#             center <- FALSE
#         }
#     } else { 
#         rm.sub <- rowMeans # a bit hacky?
#     }
#     if (center) {
#         if (verbose) { cat(" centering data by row means...") }
#         subMat <- subMat - rm.sub(subMat)  #matrix(rep(rowMeans(subMat),times=ncol(subMat)),ncol=ncol(subMat))
#         subMat[is.na(subMat)] <- 0 # replace missing with the mean
#         cat(" means for first 10 snps:\n")
#         print(round(head(rowMeans(subMat),10))) # show that centering has worked
#     } else {
#         subMat <- apply(subMat,1,row.rep)
#     }
#     if (verbose) { cat(" replaced missing data with mean (PCA cannot handle missing data)\n") }
#     #subMat <- t(subMat) # transpose
#     dimz <- dim(subMat)
#     if (!is.numeric(pcs.to.keep) | is.integer(pcs.to.keep)) { pcs.to.keep <- NA }
#     if (is.na(pcs.to.keep)) { pcs.to.keep <- min(dimz) }
#     if (pcs.to.keep > min(dimz)) { 
#         # ensure not trying to extract too many pcs
#         warning(paste0("selected too many PCs to keep [",pcs.to.keep,"], changing to ", min(dimz), "\n")) 
#         pcs.to.keep <- min(dimz)
#     } 
#     if (!SVD & (dimz[2] > dimz[1])) {
#         if (verbose) { cat(" PCA using 'princomp' (only for datasets with more samples than markers)\n") }
#         print(system.time(result <- princomp(t(subMat))))
#         PCs <- result$scores[,1:pcs.to.keep]
#         loadings <- result$loadings[,1:pcs.to.keep]
#         Evalues <- result$sdev^2 # sds are sqrt of eigenvalues
#     } else {
#         if (!SVD) {
#             if (verbose) {
#                 cat(" PCA by crossproduct and solving eigenvectors\n")
#                 cat(" obtaining crossproduct of the matrix and transpose XtX...")
#             }
#             uu <- (system.time(xtx <- crossprod(subMat)))
#             if (verbose) { 
#                 cat("took",round(uu[3]/60,1),"minutes\n")
#                 cat(" obtaining eigen vectors of the crossproduct XtX...")
#             }
#             uu <- (system.time(result <- eigen((xtx/nrow(subMat)),symmetric = TRUE)))
#             if (verbose) {  cat("took",round(uu[3]/60,1), "minutes\n") }
#             PCs <- result$vectors[,1:pcs.to.keep]
#             Evalues <- result$values
#             loadings <- NULL
#         } else {
#             pt <- "package:"; pkgset <- gsub(pt,"",search()[grep(pt,search())])
#             do.fast <- (!LAP & (all(c("irlba","bigalgebra") %in% pkgset)))
#             if (!use.bigalgebra) { do.fast <- FALSE }
#             if (verbose) {
#                 cat(" PCA by singular value decomposition...") # La.svd gives result with reversed dims. (faster?)
#             } 
#             if (return.loadings)  { nU <- pcs.to.keep } else { nU <- 0 }
#             if (!LAP) {
#                 if (do.fast) {
#                     uu <- (system.time(result <- irlba(subMat, nv = pcs.to.keep, nu = nU, mult = matmul))) 
#                 } else {
#                     if (use.bigalgebra & verbose) {
#                         warning("[without 'bigalgebra' package, PCA runs slowly for large datasets,see 'big.algebra.install.help()']\n")
#                     }
#                     uu <- (system.time(result <- svd(subMat, nv = pcs.to.keep, nu = nU)))
#                 }
#                 if (verbose) { cat("took",round(uu[3]/60,1), "minutes\n") }
#                 PCs <- result$v[,1:pcs.to.keep]
#                 #print("thus ones"); prv(result,return.loadings,nU)
#                 if (return.loadings) {
#                     loadings <- result$u[,1:pcs.to.keep]
#                 } else {
#                     loadings <- NULL
#                 }
#                 Evalues <- result$d^2 # singular values are the sqrt of eigenvalues 
#             } else {
#                 if (verbose) { cat("\n [using LAPACK alternative with La.svd]") }
#                 uu <- (system.time(result <- La.svd(subMat,nv = pcs.to.keep,nu = nU)))
#                 if (verbose) { cat("took",round(uu[3]/60,1),"minutes\n") }
#                 PCs <- t(result$vt)[,1:pcs.to.keep]  ##?
#                 # print("thOs ones")
#                 if (return.loadings) {
#                     loadings <- result$u[,1:pcs.to.keep] 
#                 } else {
#                     loadings <- NULL
#                 }
#                 Evalues <- result$d^2 # singular values are the sqrt of eigenvalues
#             }
#         }
#     }
#     rownames(PCs) <- colnames(subMat)
#     colnames(PCs) <- paste("PC",1:pcs.to.keep,sep = "")
#     if (save.pcs) {
#         ofn <- cat.path(dir$pc,pcs.fn)
#         cat(paste("~wrote PC data to file:",ofn,"\n"))
#         save(PCs,Evalues,loadings,file = ofn) }
#     if (return.loadings & exists("loadings")) {
#         colnames(loadings) <- paste("PC",1:pcs.to.keep,sep = "")
#         rownames(loadings) <- rownames(subMat)
#         out.dat <- list(PCs,Evalues,loadings)
#         names(out.dat) <- c("PCs","Evalues","loadings")
#     } else {
#         out.dat <- list(PCs,Evalues)
#         names(out.dat) <- c("PCs","Evalues")
#     }
#     return(out.dat)
# }
# 
# 
# get.big.matrix <- function(fn, dir="", verbose = FALSE)
# {
#     # loads a big.matrix either using an big.matrix description object
#     # , or this object in a binary file or text file, or points to a bigmatrix or matrix
#     if (all(dir == "")) { dir <- getwd() }
# 
#     dir.big <- dir
#     if (is.list(dir)) { if (!is.null(dir[["big"]])) { dir.big <- dir$big } }
# 
#     if (is(fn)[1] == "big.matrix.descriptor") {
#         bigMat2 <- attach.big.matrix(fn, path = dir.big)
#     } else {
#         if (is(fn)[1] == "big.matrix" | is(fn)[1] == "matrix") {
#             if (is(fn)[1] == "matrix") {
#                 bigMat2 <- as.big.matrix(fn, descriptorfile = "TEMPBIG", backingpath = dir.big)
#             } else {
#                 bigMat2 <- fn
#             }
#         } else {
#             lastchar <- substr(dir.big,nchar(dir.big),nchar(dir.big))
#             if (length(grep(".RData",fn)) == 0) {
#                 fn <- basename(fn)
#                 if (!fn %in% list.files(dir.big)) { 
#                     stop(paste("Error: big.matrix file '",fn,"' not in 'dir.big'",sep = ""))
#                 }
#                 if (verbose) { cat(" loading big matrix using text description\n") }
#                 if (lastchar == "/") { dir.big <- substr(dir.big, 1, nchar(dir.big) - 1) }
#                 bigMat2 <- attach.big.matrix(fn, path = dir.big)
#             } else {
#                 if (verbose) { cat(" loading big matrix using RData description\n") }
#                 if (lastchar != "/") { dir.big <- paste(dir.big,"/",sep = "") }
#                 filenm <- cat.path(dir.big, fn, must.exist = TRUE)
#                 dscnm <- paste(load(filenm))
#                 big.fn <- NULL
#                 for (ii in 1:length(dscnm)) {
#                     if ("big.matrix.descriptor" %in% is(get(dscnm[ii]))) {
#                         big.fn <- dscnm[ii]
#                     } 
#                 }
#                 if (!is.null(big.fn)) {
#                     descr <- get(big.fn) 
#                 } else {
#                     stop(paste("Error: didn't find bigmatrix descriptor in file",fn))
#                 }
#                 bigMat2 <- attach.big.matrix(descr, path = dir.big) 
#             }
#         }
#     }
#     return(bigMat2)
# }