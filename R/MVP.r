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


#' MVP, Memory-efficient, Visualization-enhanced, Parallel-accelerated
#'
#' Object 1: To perform GWAS using General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU model
#' Object 2: To calculate kinship among individuals using Varaden method
#' Object 3: Estimate variance components using EMMA, FaST-LMM, and HE regression
#' Object 4: Generate high-quality figures
#' 
#' @author Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param phe phenotype, n * 2 matrix, n is sample size
#' @param geno genotype, either m by n or n by m is supportable, m is marker size, n is population size
#' @param map SNP map information, SNP name, Chr, Pos
#' @param K Kinship, Covariance matrix(n * n) for random effects, must be positive semi-definite
#' @param nPC.GLM number of PCs added as fixed effects in GLM
#' @param nPC.MLM number of PCs added as fixed effects in MLM
#' @param nPC.FarmCPU number of PCs added as fixed effects in FarmCPU
#' @param CV.GLM covariates added in GLM
#' @param CV.MLM covariates added in MLM
#' @param CV.FarmCPU covariates added in FarmCPU
#' @param REML a list contains ve and vg
#' @param maxLine the number of markers handled at a time, smaller value would reduce the memory cost
#' @param ncpus number of cpus used for parallel
#' @param vc.method methods for estimating variance component("EMMA" or "HE" or "BRENT")
#' @param method the GWAS model, "GLM", "MLM", and "FarmCPU", models can be selected simutaneously, i.e. c("GLM", "MLM", "FarmCPU")
#' @param maf the threshold of minor allele frequency to filter SNPs in analysis
#' @param p.threshold if all p values generated in the first iteration are bigger than p.threshold, FarmCPU stops
#' @param QTN.threshold in second and later iterations, only SNPs with lower p-values than QTN.threshold have chances to be selected as pseudo QTNs
#' @param method.bin 'static' or 'FaST-LMM'
#' @param bin.size window size in genome
#' @param bin.selection a vector, how many windows selected
#' @param maxLoop maximum number of iterations
#' @param permutation.threshold if use a permutation cutoff or not (bonferroni cutoff)
#' @param permutation.rep number of permutation replicates
#' @param col for color of points in each chromosome on manhattan plot
#' @param memo Character. A text marker on output files
#' @param outpath the path of the output files
#' @param file.output whether to output files or not
#' @param file.type figure formats, "jpg", "tiff"
#' @param dpi resolution for output figures
#' @param threshold a cutoff line on manhattan plot, 0.05/marker size
#' @param verbose whether to print detail.
#' 
#' @export
#' @return
#' Output: MVP.return$map - SNP map information, SNP name, Chr, Pos
#' Output: MVP.return$glm.results - p-values obtained by GLM method
#' Output: MVP.return$mlm.results - p-values obtained by MLM method
#' Output: MVP.return$farmcpu.results - p-values obtained by FarmCPU method
#'
#' @examples
#' \donttest{
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' mapPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.map", package = "rMVP")
#' map <- read.table(mapPath , head = TRUE)
#' 
#' opts <- options(rMVP.OutputLog2File = FALSE)
#' 
#' mvp <- MVP(phe=phenotype, geno=genotype, map=map, maxLoop=3,
#'   method=c("GLM", "MLM", "FarmCPU"), file.output=FALSE, ncpus=1)
#' str(mvp)
#' 
#' options(opts)
#' }
MVP <-
function(phe, geno, map, K=NULL, nPC.GLM=NULL, nPC.MLM=NULL, nPC.FarmCPU=NULL,
         CV.GLM=NULL, CV.MLM=NULL, CV.FarmCPU=NULL, REML=NULL, maxLine=10000, 
         ncpus=detectCores(logical = FALSE), vc.method=c("BRENT", "EMMA", "HE"), 
         method=c("GLM", "MLM", "FarmCPU"), maf=NULL, p.threshold=NA, 
         QTN.threshold=0.01, method.bin="static", bin.size=c(5e5,5e6,5e7), 
         bin.selection=seq(10,100,10), maxLoop=10, permutation.threshold=FALSE, 
         permutation.rep=100, memo=NULL, outpath=getwd(),
         col=c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", 
         "#d581b7", "#83d3ad", "#7c162c", "#26755d"), file.output=TRUE, 
         file.type="jpg", dpi=300, threshold=0.05, verbose=TRUE
) {

    # Compatible with old ways
    if (is.logical(file.output)) {
        if (file.output == TRUE) {
            file.output <- c("pmap", "pmap.signal", "plot", "log")
        } else if (file.output == FALSE) {
            file.output <- c()
        }
    }

    for(mt in method){
        if(!mt %in% c("GLM", "MLM", "FarmCPU"))
            stop("Unknow method: ", mt)
    }
    
    # Set output path of log file 
    logging.outpath <- NULL
    if ("log" %in% file.output) {
      logging.outpath <- outpath
    }
    logging.initialize("MVP", logging.outpath)

    MVP.Version(width = 65, verbose = verbose)
    time_start <- Sys.time()
    logging.log("Start:", format(time_start, format = "%F %T %Z"), "\n", verbose = verbose)
    if ("log" %in% file.output) {
        logging.log("The log has been output to the file:", get("logging.file", envir = package.env), "\n", verbose = verbose)
    }
    vc.method <- match.arg(vc.method)
    if (nrow(phe) != ncol(geno) & nrow(phe) != nrow(geno)) stop("The number of individuals in phenotype and genotype doesn't match!")
    if (nrow(map) != ncol(geno) & nrow(map) != nrow(geno)) stop("The number of markers in genotype and map doesn't match!")
    if (!is.big.matrix(geno))    stop("genotype should be in 'big.matrix' format.")
    
    #list -> matrix
    map <- as.data.frame(map)
    for(i in 1 : ncol(map)){
        if(is.factor(map[, i])) map[, i] <- as.character.factor(map[, i])
    }

    na.index <- NULL
    if (!is.null(CV.GLM)) {
        CV.GLM <- as.matrix(CV.GLM)
    if (nrow(CV.GLM) != nrow(phe)) stop("The number of individuals in covariates and phenotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.GLM), arr.ind = TRUE)[, 1])
    }
    if (!is.null(CV.MLM)) {
        CV.MLM <- as.matrix(CV.MLM)
        if (nrow(CV.MLM) != nrow(phe)) stop("The number of individuals in covariates and phenotype doesn't match!")
       na.index <- c(na.index, which(is.na(CV.MLM), arr.ind = TRUE)[, 1])
    }
    if (!is.null(CV.FarmCPU)) {
        CV.FarmCPU <- as.matrix(CV.FarmCPU)
        if (nrow(CV.FarmCPU) != nrow(phe)) stop("The number of individuals in covariates and phenotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.FarmCPU), arr.ind = TRUE)[, 1])
    }
    na.index <- unique(na.index)

    #Data information
    MrkByCol <- nrow(phe) == nrow(geno)
    m <- ifelse(MrkByCol, ncol(geno), nrow(geno))
    n <- nrow(phe)
    logging.log(paste("Input data has", n, "individuals and", m, "markers"), "\n", verbose = verbose)
    logging.log(paste("Markers are detected to be stored by", ifelse(MrkByCol, "column", "row")), "\n", verbose = verbose)
    logging.log("Analyzed trait:", colnames(phe)[2], "\n", verbose = verbose)
    logging.log("Number of threads used:", ncpus, "\n", verbose = verbose)
    hpclib <- grepl("mkl", sessionInfo()$LAPACK) | grepl("openblas", sessionInfo()$LAPACK) | eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))
	if(!hpclib){
		logging.log("No high performance math library detected! The computational efficiency would be greatly reduced\n", verbose = verbose)
	}else{
		if(grepl("mkl", sessionInfo()$LAPACK) | eval(parse(text = "!inherits(try(Revo.version,silent=TRUE),'try-error')"))){
			logging.log("Math Kernel Library is detected, nice job!\n", verbose = verbose)
		}else{
			logging.log("OpenBLAS Library is detected, nice job!\n", verbose = verbose)
		}
	}

    #remove samples with missing phenotype
    seqTaxa <- which(!is.na(phe[,2]))
    if (length(na.index) != 0) seqTaxa <- intersect(seqTaxa, c(1:n)[-na.index])
    if (length(seqTaxa) == 0)   stop("no effective individuals left due to missings")
    if (length(seqTaxa) == n) seqTaxa <- NULL
    if (!is.null(seqTaxa)) {
        logging.log("Total", n - length(seqTaxa), "individuals are removed due to missings", "\n", verbose = verbose)
        phe = phe[seqTaxa,]
        if (!is.null(K)) { K = K[seqTaxa, seqTaxa] }
        if (!is.null(CV.GLM)) { CV.GLM = CV.GLM[seqTaxa, , drop = FALSE] }
        if (!is.null(CV.MLM)) { CV.MLM = CV.MLM[seqTaxa, , drop = FALSE] }
        if (!is.null(CV.FarmCPU)) { CV.FarmCPU = CV.FarmCPU[seqTaxa, , drop = FALSE] }
        if(length(seqTaxa) < n * 0.8){
            logging.log("Re-build memory-mapping file for remaining individuals", "\n", verbose = verbose)
            if(!MrkByCol){
                geno <- deepcopy(geno, cols = seqTaxa)
            }else{
                geno <- deepcopy(geno, rows = seqTaxa)
            }
            seqTaxa <- NULL
        }
    }

    logging.log("Calculate allele frequency...", "\n", verbose = verbose)
    marker_freq <- BigRowMean(geno@address, MrkByCol, threads = ncpus, geno_ind = seqTaxa) / 2
    map$MAF <- ifelse(marker_freq > 0.5, 1 - marker_freq, marker_freq)

    #remove SNPs with low MAF
    geno_marker_index <- NULL
    map_sub <- map
    if (!is.null(maf)){
        if(length(maf) != 1) stop("maf should be a value")
        if(maf <= 0 || maf >= 0.5) stop("maf should be at the range of 0-0.5")
        geno_marker_index <- which(map$MAF >= maf)
        if(length(geno_marker_index) == 0) stop(paste("MAFs of all markers are smaller than the threshold", maf))
        if(length(geno_marker_index) == 1) stop(paste("only 1 marker left on the given MAF threshold", maf))
        if(length(geno_marker_index) == m)  geno_marker_index <- NULL
    }
    if(!is.null(geno_marker_index)){
        logging.log("Total", m - length(geno_marker_index), "markers are removed at MAF threshold", maf, "\n", verbose = verbose)
        m <- length(geno_marker_index)
        map_sub <- map[geno_marker_index, ]
        marker_freq <- marker_freq[geno_marker_index]
        if(length(geno_marker_index) < m * 0.8){
            if(!is.null(seqTaxa)){
                logging.log("Re-build memory-mapping file for remaining individuals and markers", "\n", verbose = verbose)
                if(MrkByCol){
                    geno <- deepcopy(geno, rows = seqTaxa, cols = geno_marker_index)
                }else{
                    geno <- deepcopy(geno, cols = seqTaxa, rows = geno_marker_index)
                }
                seqTaxa <- NULL
            }else{
                logging.log("Re-build memory-mapping file for remaining markers", "\n", verbose = verbose)
                if(MrkByCol){
                    geno <- deepcopy(geno, cols = geno_marker_index)
                }else{
                    geno <- deepcopy(geno, rows = geno_marker_index)
                }
            }
            geno_marker_index <- NULL
        }
    }

    #initial results
    glm.results <- NULL
    mlm.results <- NULL
    farmcpu.results <- NULL
    
    #indicators for models
    glm.run <- "GLM" %in% method
    mlm.run <- "MLM" %in% method
    farmcpu.run <- "FarmCPU" %in% method
    
    nPC <- suppressWarnings(max(nPC.GLM, nPC.MLM, nPC.FarmCPU, na.rm = TRUE))
    if (nPC <= 0) {
        nPC <- NULL
    } else if (nPC < 3) {
        nPC <- 3
    }

    if (!is.null(K)) { K <- as.matrix(K) }
    if (!is.null(nPC) | "MLM" %in% method) {
        if (is.null(K)) {
            K <- MVP.K.VanRaden(
                M = geno, 
                ind_idx = seqTaxa, 
                mrk_idx = geno_marker_index,
                mrk_freq = marker_freq,
                mrk_bycol = MrkByCol,
                maxLine = maxLine, 
                cpu = ncpus, 
              verbose = verbose,
              checkNA = FALSE
            )
        }
        logging.log("Eigen Decomposition on GRM", "\n", verbose = verbose)
        eigenK <- eigen(K, symmetric = TRUE)
        if (!is.null(nPC)) {
            ipca <- eigenK$vectors[, 1:nPC]
            logging.log("Deriving PCs successfully", "\n", verbose = verbose)
        }
        if (("MLM" %in% method) & vc.method == "BRENT") { K <- NULL; gc()}
        if (!"MLM" %in% method) { rm(eigenK); rm(K); gc() }
    }

    if (!is.null(nPC)) {
        #CV for GLM
        if (glm.run) {
            if (!is.null(CV.GLM)) {
                logging.log("Number of provided covariates of GLM:", ncol(CV.GLM), "\n", verbose = verbose)
                if (!is.null(nPC.GLM)) {
                    logging.log("Number of PCs included:", nPC.GLM, "\n", verbose = verbose)
                    CV.GLM <- cbind(ipca[,1:nPC.GLM], CV.GLM)
                }
            } else if (!is.null(nPC.GLM)) {
                logging.log("Number of PCs included in GLM:", nPC.GLM, "\n", verbose = verbose)
                CV.GLM <- ipca[,1:nPC.GLM, drop = FALSE]
            }
        }
        
        #CV for MLM
        if (mlm.run) {
            if (!is.null(CV.MLM)) {
                logging.log("Number of provided covariates of MLM:", ncol(CV.MLM), "\n", verbose = verbose)
                if (!is.null(nPC.MLM)) {
                    logging.log("Number of PCs included:", nPC.MLM, "\n", verbose = verbose)
                    CV.MLM <- cbind(ipca[,1:nPC.MLM], CV.MLM)
                }
            } else if (!is.null(nPC.MLM)) {
                logging.log("Number of PCs included in MLM:", nPC.MLM, "\n", verbose = verbose)
                CV.MLM <- ipca[,1:nPC.MLM, drop = FALSE]
            }
        }
        
        #CV for FarmCPU
        if (farmcpu.run) {
            if (!is.null(CV.FarmCPU)) {
                logging.log("Number of provided covariates of FarmCPU:", ncol(CV.FarmCPU), "\n", verbose = verbose)
                if (!is.null(nPC.FarmCPU)) {
                    logging.log("Number of PCs included:", nPC.FarmCPU, "\n", verbose = verbose)
                    CV.FarmCPU <- cbind(ipca[,1:nPC.FarmCPU], CV.FarmCPU)
                }
            }else if (!is.null(nPC.FarmCPU)) { 
                logging.log("Number of PCs included in FarmCPU:", nPC.FarmCPU, "\n", verbose = verbose)
                CV.FarmCPU <- ipca[,1:nPC.FarmCPU]
            }
        }  
    }else{
        if (glm.run) {
            if (!is.null(CV.GLM)) {
                logging.log("Number of provided covariates of GLM:", ncol(CV.GLM), "\n", verbose = verbose)
            }
        }
        if (mlm.run) {
            if (!is.null(CV.MLM)) {
                logging.log("Number of provided covariates of MLM:", ncol(CV.MLM), "\n", verbose = verbose)
            }
        }
        if (farmcpu.run) {
            if (!is.null(CV.FarmCPU)) {
                logging.log("Number of provided covariates of FarmCPU:", ncol(CV.FarmCPU), "\n", verbose = verbose)
            }
        }
    }

    #GWAS
    logging.log("-------------------------GWAS Start-------------------------", "\n", verbose = verbose)
    if (glm.run) {
        logging.log("General Linear Model (GLM) Start...", "\n", verbose = verbose)
        glm.results <- MVP.GLM(phe=phe, geno=geno, CV=CV.GLM, ind_idx=seqTaxa, mrk_idx=geno_marker_index, mrk_bycol = MrkByCol, maxLine = maxLine, cpu=ncpus, verbose = verbose);gc()
        colnames(glm.results) <- c("Effect", "SE", paste(colnames(phe)[2],"GLM",sep="."))
        z = glm.results[, 1]/glm.results[, 2]
        lambda = median(z^2, na.rm=TRUE)/qchisq(1/2, df = 1,lower.tail=FALSE)
        logging.log("Genomic inflation factor (lambda):", round(lambda, 4), "\n", verbose = verbose)
        if ("pmap" %in% file.output) {
            logging.log("Writing results to local file", "\n", verbose = verbose)
            write.csv(x = cbind(map_sub, glm.results), 
                    file = file.path(outpath, paste(colnames(phe)[2], ".GLM.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                    row.names = FALSE)
        }
    }

    if (mlm.run) {
        logging.log("Mixed Linear Model (MLM) Start...", "\n", verbose = verbose)
        mlm.results <- MVP.MLM(phe=phe, geno=geno, K=K, eigenK=eigenK, CV=CV.MLM, ind_idx=seqTaxa, mrk_idx=geno_marker_index, mrk_bycol = MrkByCol, maxLine = maxLine, cpu=ncpus, vc.method=vc.method, verbose = verbose);gc()
        colnames(mlm.results) <- c("Effect", "SE", paste(colnames(phe)[2],"MLM",sep="."))
        z = mlm.results[, 1]/mlm.results[, 2]
        lambda = median(z^2, na.rm=TRUE)/qchisq(1/2, df = 1,lower.tail=FALSE)
        logging.log("Genomic inflation factor (lambda):", round(lambda, 4), "\n", verbose = verbose)
        if ("pmap" %in% file.output) {
            logging.log("Writing results to local file", "\n", verbose = verbose)
            write.csv(x = cbind(map_sub, mlm.results), 
                    file = file.path(outpath, paste(colnames(phe)[2], ".MLM.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                    row.names = FALSE)
        }
    }
    
    if (farmcpu.run) {
        logging.log("FarmCPU Start...", "\n", verbose = verbose)
        farmcpu.results <- MVP.FarmCPU(phe=phe, geno=geno, map=map[,1:3], CV=CV.FarmCPU, ind_idx=seqTaxa, mrk_idx=geno_marker_index, maxLine = maxLine, ncpus=ncpus, memo="MVP.FarmCPU", p.threshold=p.threshold, QTN.threshold=QTN.threshold, method.bin=method.bin, bin.size=bin.size, bin.selection=bin.selection, maxLoop=maxLoop, verbose = verbose)
        colnames(farmcpu.results) <- c("Effect", "SE", paste(colnames(phe)[2],"FarmCPU",sep="."))
        z = farmcpu.results[, 1]/farmcpu.results[, 2]
        lambda = median(z^2, na.rm=TRUE)/qchisq(1/2, df = 1,lower.tail=FALSE)
        logging.log("Genomic inflation factor (lambda):", round(lambda, 4), "\n", verbose = verbose)
        if ("pmap" %in% file.output) {
            logging.log("Writing results to local file", "\n", verbose = verbose)
            write.csv(x = cbind(map_sub, farmcpu.results), 
                    file = file.path(outpath, paste(colnames(phe)[2], ".FarmCPU.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                    row.names = FALSE)
        }
    }
    
    MVP.return <- list(map=map_sub, glm.results=glm.results, mlm.results=mlm.results, farmcpu.results=farmcpu.results)
    
    if(permutation.threshold){
        # set.seed(12345)
        i=1
        for(i in 1:permutation.rep){
            index = 1:nrow(phe)
            index.shuffle = sample(index,length(index),replace=FALSE)
            myY.shuffle = phe
            myY.shuffle[,2] = myY.shuffle[index.shuffle,2]
            #GWAS using t.test...
            myPermutation = MVP.GLM(phe=myY.shuffle[,c(1,2)], geno=geno, ind_idx=seqTaxa, mrk_idx=geno_marker_index, maxLine=maxLine, cpu=ncpus)
            pvalue = min(myPermutation[,3],na.rm=TRUE)
            if(i==1){
                pvalue.final=pvalue
            }else{
                pvalue.final=c(pvalue.final,pvalue)
            }
        }#end of permutation.rep
        permutation.cutoff = sort(pvalue.final)[ceiling(permutation.rep*0.05)]
        threshold = permutation.cutoff * m
    }
    logging.log(paste0("Significant level: ", formatC(threshold/m, format = "e", digits = 2)), "\n", verbose = verbose)
    if ("pmap.signal" %in% file.output) {
        if (glm.run) {
            index <- which(glm.results[, ncol(glm.results)] < threshold/m)
            if (length(index) != 0) {
              write.csv(x = cbind.data.frame(map_sub, glm.results)[index, ], 
                        file = file.path(outpath, paste(colnames(phe)[2], ".GLM_signals.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                        row.names = FALSE)
            }
        }
        if (mlm.run) {
            index <- which(mlm.results[, ncol(mlm.results)] < threshold/m)
            if (length(index) != 0) {
              write.csv(x = cbind.data.frame(map_sub, mlm.results)[index, ], 
                        file = file.path(outpath, paste(colnames(phe)[2], ".MLM_signals.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                        row.names = FALSE)
            }
        }
        if (farmcpu.run) {
            index <- which(farmcpu.results[, ncol(farmcpu.results)] < threshold/m)
            if (length(index) != 0) {
              write.csv(x = cbind.data.frame(map_sub, farmcpu.results)[index, ], 
                        file = file.path(outpath, paste(colnames(phe)[2], ".FarmCPU_signals.", memo, ifelse(is.null(memo),"csv",".csv"), sep = "")),
                        row.names = FALSE)
            }
        }
    }
    if ("plot" %in% file.output) {
        logging.log("---------------------Visualization Start--------------------", "\n", verbose = verbose)
        logging.log("Phenotype distribution Plotting", "\n", verbose = verbose)
        MVP.Hist(memo=memo, outpath=outpath, file.output=TRUE, phe=phe, file.type=file.type, col=col, dpi=dpi)
        #plot3D <- !is(try(library("rgl"),silent=TRUE), "try-error")
        plot3D <- FALSE
        if(!is.null(nPC)){
          MVP.PCAplot(
            ipca[,1:3],
            col=col,
            plot3D=plot3D,
            file.output=TRUE,
            file.type=file.type,
            outpath=outpath, 
            memo = ifelse(is.null(memo), colnames(phe)[2], paste(colnames(phe)[2], memo, sep=".")),
            dpi=dpi,
          )
        }
        
        MVP.Report(
            MVP.return,
            col=col,
            plot.type=c("c","m","q","d"),
            file.output=TRUE,
            file.type=file.type,
            outpath=outpath, 
            memo = memo,
            chr.den.col=c("darkgreen", "yellow", "red"),
            dpi=dpi,
            threshold=threshold/m,
        )

        if(sum(c(is.null(glm.results), is.null(mlm.results), is.null(farmcpu.results))) < 2) {
            MVP.Report(
                MVP.return,
                col=col,
                plot.type=c("m","q"),
                multracks=TRUE,
                file.output=TRUE,
                file.type=file.type,
                outpath=outpath, 
                memo = memo,
                dpi=dpi,
                threshold=threshold/m
            )
        }
    }
    time_end <- Sys.time()
    if (length(file.output) > 0) {
      logging.log("Results are stored at Working Directory:", outpath, "\n", verbose = verbose)
    }
    logging.log("End:", format(time_end, format = "%F %T %Z"), "\n", verbose = verbose)
    time_diff <- as.numeric(time_end) - as.numeric(time_start)
    h <- time_diff %/% 3600
	m <- (time_diff %% 3600) %/% 60
	s <- ((time_diff %% 3600) %% 60)
	index <- which(c(h, m, s) != 0)
	num <- c(h, m, s)[index]
	num <- round(num, 0)
	char <- c("h", "m", "s")[index]
	logging.log("Total running time:", paste(num, char, sep="", collapse=""), "\n", verbose = verbose)
    print_accomplished(width = 60, verbose = verbose)
    
    return(invisible(MVP.return))
}#end of MVP function
