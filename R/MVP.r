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


#' MVP, Memory-efficient, Visualization-enhanced, Parallel-accelerated
#'
#' Object 1: To perform GWAS using General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU model
#' Object 2: To calculate kinship among individuals using Varaden method
#' Object 3: Estimate variance components using EMMA, FaST-LMM, and HE regression
#' Object 4: Generate high-quality figures
#' 
#' Build date: Aug 30, 2017
#' Last update: Dec 14, 2018
#' 
#' @author Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param phe phenotype, n * 2 matrix, n is sample size
#' @param geno Genotype in bigmatrix format; m * n, m is marker size, n is sample size
#' @param map SNP map information, SNP name, Chr, Pos
#' @param K Kinship, Covariance matrix(n * n) for random effects, must be positive semi-definite
#' @param nPC.GLM number of PCs added as fixed effects in GLM
#' @param nPC.MLM number of PCs added as fixed effects in MLM
#' @param nPC.FarmCPU number of PCs added as fixed effects in FarmCPU
#' @param perc percentage of total SNPs selected for PCA
#' @param CV.GLM covariates added in GLM
#' @param CV.MLM covariates added in MLM
#' @param CV.FarmCPU covariates added in FarmCPU
#' @param REML a list contains ve and vg
#' @param priority speed or memory
#' @param ncpus number of cpus used for parallel
#' @param vc.method methods for estimating variance component("EMMA" or "GEMMA")
#' @param method the GWAS model, "GLM", "MLM", and "FarmCPU", models can be selected simutaneously, i.e. c("GLM", "MLM", "FarmCPU")
#' @param maxLine when the priority is 'memory', users can change this parameter to limit the memory
#' @param memo a marker added on output file name
#' @param P a start p value for each SNP
#' @param method.sub, method.sub.final method used in substitution process
#' @param method.sub.final method used in substitution process, five options: 'penalty', 'reward', 'mean', 'median', or 'onsite'
#' @param method.bin EMMA or FaSTLMM
#' @param bin.size window size in genome
#' @param bin.selection a vector, how many windows selected
#' @param Prior four columns, SNP name, Chr, Pos, P
#' @param maxLoop maximum number of iterations
#' @param threshold.output output GWAS results only for SNPs with p value lower than the threshold.output
#' @param p.threshold if all p values in the 1st iteration are bigger than p.threshold, FarmCPU stops
#' @param QTN.threshold Only SNPs have a more significant p value than QTN.threshold have chance to be selected as pseudo QTNs
#' @param bound maximum number of SNPs selected as pseudo QTNs for each iteration
#' @param outward the direction of circular Manhattan plot
#' @param permutation.threshold if use a permutation cutoff or not (bonferroni cutoff)
#' @param permutation.rep number of permutation replicates
#' @param bar if TRUE, the progress bar will be drawn on the terminal
#' @param col for color of points in each chromosome on manhattan plot
#' @param plot.type "b" (both Manhattan plot and qq plot will be draw) or "q" (qq plot only)
#' @param file.output whether to output files or not
#' @param file figure formats, "jpg", "tiff"
#' @param dpi resolution for output figures
#' @param threshold a cutoff line on manhattan plot, 0.05/marker size
#' @param Ncluster number of colors used for drawing PC 1 and PC 2
#' @param signal.cex point size on output figures
#' @param box logical, if TRUE, the box frame will be added in output figure
#'
#' @export
#' @return a m * 2 matrix, the first column is the SNP effect, the second column is the P values
#' Output: MVP.return$map - SNP map information, SNP name, Chr, Pos
#' Output: MVP.return$glm.results - p-values obtained by GLM method
#' Output: MVP.return$mlm.results - p-values obtained by MLM method
#' Output: MVP.return$farmcpu.results - p-values obtained by FarmCPU method
#'
#' @examples
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header=TRUE)
#' print(dim(phenotype))
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' print(dim(genotype))
#' mapPath <- system.file("extdata", "07_other", "mvp.map", package = "rMVP")
#' map <- read.table(mapPath , head = TRUE)
#' mvp <- MVP(phe=phenotype, geno=genotype, map=map, 
#'   method=c("GLM", "MLM", "FarmCPU"), file.output=FALSE, ncpus=1)
#' str(mvp)
MVP <-
function(phe, geno, map, K=NULL, nPC.GLM=NULL, nPC.MLM=NULL, nPC.FarmCPU=NULL, perc=1, CV.GLM=NULL, CV.MLM=NULL, CV.FarmCPU=NULL, REML=NULL, priority="speed", ncpus=detectCores(logical = FALSE), vc.method="EMMA", method="MLM", maxLine=1000, memo=NULL, P=NULL, method.sub="reward", method.sub.final="reward", method.bin="static", bin.size=c(5e5,5e6,5e7), bin.selection=seq(10,100,10), Prior=NULL, maxLoop=10, threshold.output=1, p.threshold=NA, QTN.threshold=0.01, bound=NULL, outward=FALSE,
permutation.threshold=FALSE, permutation.rep=100, bar=TRUE, col=c("dodgerblue4","olivedrab4","violetred","darkgoldenrod1","purple4"), plot.type="b", file.output=TRUE, file="jpg", dpi=300, threshold=0.05, Ncluster=1, signal.cex=0.8, box=FALSE
) {
    if (Sys.info()[['sysname']] == 'Windows')
        ncpus <- 1
    
    MVP.Version(width = 60)
    if(nrow(phe) != ncol(geno)) stop("The number of individuals in phenotype and genotype doesn't match!")
    #list -> matrix
    map <- as.matrix(map)
    na.index <- NULL
    if(!is.null(CV.GLM)){
        CV.GLM <- as.matrix(CV.GLM)
    if(nrow(CV.GLM) != ncol(geno)) stop("The number of individuals in covariates and genotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.GLM), arr.ind=TRUE)[, 1])
    }
    if(!is.null(CV.MLM)){
        CV.MLM <- as.matrix(CV.MLM)
        if(nrow(CV.MLM) != ncol(geno)) stop("The number of individuals in covariates and genotype doesn't match!")
       na.index <- c(na.index, which(is.na(CV.MLM), arr.ind=TRUE)[, 1])
    }
    if(!is.null(CV.FarmCPU)){
        CV.FarmCPU <- as.matrix(CV.FarmCPU)
        if(nrow(CV.FarmCPU) != ncol(geno)) stop("The number of individuals in covariates and genotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.FarmCPU), arr.ind=TRUE)[, 1])
    }
    na.index <- unique(na.index)
    
    #remove samples with missing phenotype
    seqTaxa = which(!is.na(phe[,2]))
    if(length(na.index) != 0) seqTaxa <- intersect(seqTaxa, seq_len(nrow(phe))[-na.index])
    #file.exsits()
    if(length(seqTaxa) != length(phe[,2])){
        try(unlink(c("geno.temp.bin","geno.temp.desc")), silent=TRUE)
        geno = deepcopy(geno, cols=seqTaxa, backingfile="geno.temp.bin", descriptorfile="geno.temp.desc")
        phe = phe[seqTaxa,]
        if(!is.null(K)){K = K[seqTaxa, seqTaxa]}
        if(!is.null(CV.GLM)){CV.GLM = CV.GLM[seqTaxa,]}
        if(!is.null(CV.MLM)){CV.MLM = CV.MLM[seqTaxa,]}
        if(!is.null(CV.FarmCPU)){CV.FarmCPU = CV.FarmCPU[seqTaxa,]}
        rm(geno)
            gc()
            genoName <- "geno.temp.desc"
        geno <- attach.big.matrix(genoName)
        unlink(c("geno.temp.desc", "geno.temp.bin"))
    }
    #Data information
    m=nrow(geno)
    n=ncol(geno)
    print(paste("Input data has", n, "individuals,", m, "markers", sep=" "))
    
    #initial results
    glm.results <- NULL
    mlm.results <- NULL
    farmcpu.results <- NULL
    
    #indicators for models
    glm.run <- "GLM" %in% method
    mlm.run <- "MLM" %in% method
    farmcpu.run <- "FarmCPU" %in% method
    
    if(!is.null(nPC.GLM)|!is.null(nPC.MLM)|!is.null(nPC.FarmCPU)){
        nPC <- max(nPC.GLM, nPC.MLM, nPC.FarmCPU)
        #if(nPC < 3){
        #nPC <- 3
        #}
    }else{
        nPC <- NULL
    }
    
    if (!is.null(nPC)) {
        ipca <- MVP.PCA(M = geno, perc = perc, pcs.keep = nPC)$PCs
        if (file.output) {
            filebck <- paste0("MVP.", colnames(phe)[2], memo, ".pc.bin")
            filedes <- paste0("MVP.", colnames(phe)[2], memo, ".pc.desc")
            
            if (file.exists(filebck)) file.remove(filebck)
            if (file.exists(filedes)) file.remove(filedes)
            PC.backed <- big.matrix(
                nrow = nrow(ipca), 
                ncol = ncol(ipca), 
                type = "double", 
                backingfile = filebck,
                descriptorfile = filedes
            )
            
            PC.backed[, ] <- ipca[, ]
            flush(PC.backed)
            rm(list = c("PC.backed"))
        }
        #CV for GLM
        if(glm.run){
            if(!is.null(CV.GLM)){
                if(!is.null(nPC.GLM)){
                    CV.GLM <- cbind(ipca[,seq_len(nPC.GLM)], CV.GLM)
                }
            }else{
                CV.GLM <- ipca[,seq_len(nPC.GLM)]
            }
        }
        
        #CV for MLM
        if(mlm.run){
            if(!is.null(CV.MLM)){
                if(!is.null(nPC.MLM)){
                    CV.MLM <- cbind(ipca[,seq_len(nPC.MLM)], CV.MLM)
                }
            }else{
                CV.MLM <- ipca[,seq_len(nPC.MLM)]
            }
        }
        
        #CV for FarmCPU
        if(farmcpu.run){
            if(!is.null(CV.FarmCPU)){
                if(!is.null(nPC.FarmCPU)){
                    CV.FarmCPU <- cbind(ipca[,seq_len(nPC.FarmCPU)], CV.FarmCPU)
                }
            }else{
                CV.FarmCPU <- ipca[,seq_len(nPC.FarmCPU)]
            }
        }
        
    }
  
    #GWAS
    print("GWAS Start...")
    # ------------------------- GLM -------------------------
    if (glm.run) {
        print("General Linear Model (GLM) Start ...")

        glm.results <-
            MVP.GLM(
                phe = phe,
                geno = geno,
                CV = CV.GLM,
                cpu = ncpus,
                memo = "MVP.GLM",
                bar = bar
            )
        gc()
        colnames(glm.results) <-
            c("effect", paste(colnames(phe)[2], "GLM", sep = "."))
        if (file.output)
            write.csv(
                cbind(map, glm.results),
                paste("MVP.", colnames(phe)[2], ".GLM", ".csv", sep = ""),
                row.names = FALSE
            )
    }
    
    # ------------------------- MLM -------------------------
    if (mlm.run) {
        print("Mixed Linear Model (MLM) Start ...")
        mlm.results <-
            MVP.MLM(
                phe = phe,
                geno = geno,
                K = K,
                priority = priority,
                CV = CV.MLM,
                cpu = ncpus,
                bar = bar,
                maxLine = maxLine,
                vc.method = vc.method,
                file.output = file.output,
                memo = "MVP.MLM"
            )
        gc()
        colnames(mlm.results) <-
            c("effect", paste(colnames(phe)[2], "MLM", sep = "."))
        if (file.output)
            write.csv(
                cbind(map, mlm.results),
                paste("MVP.", colnames(phe)[2], ".MLM", ".csv", sep = ""),
                row.names = FALSE
            )
    }
    
    # ----------------------- FarmCPU -----------------------
    if (farmcpu.run) {
        print("FarmCPU Start ...")
        farmcpu.results <-
            MVP.FarmCPU(
                phe = phe,
                geno = geno,
                map = map,
                priority = priority,
                CV = CV.FarmCPU,
                ncpus = ncpus,
                bar = bar,
                memo = "MVP.FarmCPU",
                P = P,
                method.sub = method.sub,
                method.sub.final = method.sub.final,
                method.bin = method.bin,
                bin.size = bin.size,
                bin.selection = bin.selection,
                Prior = Prior,
                maxLoop = maxLoop,
                threshold.output = threshold.output,
                p.threshold = p.threshold,
                QTN.threshold = QTN.threshold
            )
        colnames(farmcpu.results) <-
            c("effect", paste(colnames(phe)[2], "FarmCPU", sep = "."))
        if (file.output)
            write.csv(
                cbind(map, farmcpu.results),
                paste("MVP.", colnames(phe)[2], ".FarmCPU", ".csv", sep = ""),
                row.names = FALSE
            )
    }
    
    MVP.return <-
        list(
            map = map,
            glm.results = glm.results,
            mlm.results = mlm.results,
            farmcpu.results = farmcpu.results
        )
    
    if (permutation.threshold) {
        pvalue.final <- c()
        for (i in seq_len(permutation.rep)) {
            index         <- seq_len(nrow(phe))
            index.shuffle <- sample(index, length(index), replace = FALSE)
            myY.shuffle   <- phe
            myY.shuffle[, 2] <- myY.shuffle[index.shuffle, 2]
            
            #GWAS using t.test...
            myPermutation <-
                MVP.GLM(phe = myY.shuffle[, c(1, 2)], geno = geno, cpu = ncpus)
            
            pvalue       <- min(myPermutation[, 2], na.rm = TRUE)
            pvalue.final <- c(pvalue.final,pvalue)
        }#end of permutation.rep
        permutation.cutoff <- sort(pvalue.final)[ceiling(permutation.rep*0.05)]
        threshold <- permutation.cutoff * m
    }
    
    print(paste("Significance Level: ", threshold/m, sep=""))
    if(file.output){
        print("Visualization Start...")
        print("Phenotype distribution Plotting...")
        MVP.Hist(phe=phe, file.type=file, col=col, dpi=dpi)
        #plot3D <- !is(try(library("rgl"),silent=TRUE), "try-error")
        plot3D <- TRUE
        if(!is.null(nPC)){
            MVP.PCAplot(
                ipca[,1:3],
                col=col,
                plot3D=plot3D,
                Ncluster=Ncluster,
                file=file,
                dpi=dpi,
                box=box
            )
        }
        
        MVP.Report(
            MVP.return,
            col=col,
            box=box,
            plot.type=c("c","m","q","d"),
            file.output=TRUE,
            file=file,
            dpi=dpi,
            threshold=threshold/m,
            signal.cex=signal.cex,
            outward=outward
        )
        
        if(sum(c(is.null(glm.results), is.null(mlm.results), is.null(farmcpu.results))) < 2) {
            MVP.Report(
                MVP.return,
                col=col,
                plot.type=c("m","q"),
                multracks=TRUE,
                outward=outward,
                file.output=TRUE,
                file=file,
                dpi=dpi,
                box=box,
                threshold=threshold/m,
                signal.cex=signal.cex
            )
        }
    }
    print_accomplished(width = 60)
    return(MVP.return)
}#end of MVP function
