MVP <-
function(phe, geno, map, K=NULL, nPC.GLM=NULL, nPC.MLM=NULL, nPC.FarmCPU=NULL, perc=1, CV.GLM=NULL, CV.MLM=NULL, CV.FarmCPU=NULL, REML=NULL, priority="speed", ncpus=detectCores(logical = F), vc.method="EMMA", method="MLM", maxLine=1000, memo=NULL, P=NULL, method.sub="reward", method.sub.final="reward", method.bin="static", bin.size=c(5e5,5e6,5e7), bin.selection=seq(10,100,10), Prior=NULL, maxLoop=10, threshold.output=1, iteration.output=FALSE, p.threshold=NA, QTN.threshold=NULL, bound=NULL, outward=FALSE,
permutation.threshold=FALSE, permutation.rep=100, bar=TRUE, col=c("dodgerblue4","olivedrab4","violetred","darkgoldenrod1","purple4"), plot.type="b", file.output=TRUE, file="jpg", dpi=300, threshold=0.05, Ncluster=1, signal.cex=0.8, box=FALSE
)
{
###################################################################################################
# Object:  MVP, Memory-efficient, Visualization-enhanced, Parallel-accelerated
#
# Object 1: To perform GWAS using General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU model
# Object 2: To calculate kinship among individuals using Varaden method
# Object 3: Estimate variance components using EMMA, FaST-LMM, and HE regression
# Object 4: Generate high-quality figures
#
# Input:
# phe: phenotype, n * 2 matrix, n is sample size
# geno: Genotype in bigmatrix format; m * n, m is marker size, n is sample size
# map: SNP map information, SNP name, Chr, Pos
# K: Kinship, Covariance matrix(n * n) for random effects, must be positive semi-definite
# nPC.GLM: number of PCs added as fixed effects in GLM
# nPC.MLM: number of PCs added as fixed effects in MLM
# nPC.FarmCPU: number of PCs added as fixed effects in FarmCPU
# perc: percentage of total SNPs selected for PCA
# CV.GLM: covariates added in GLM
# CV.MLM: covariates added in MLM
# CV.FarmCPU: covariates added in FarmCPU
# REML: a list contains ve and vg
# priority: speed or memory
# ncpuc: number of cpus used for parallel
# vc.method: methods for estimating variance component("EMMA" or "GEMMA")
# method: the GWAS model, "GLM", "MLM", and "FarmCPU", models can be selected simutaneously, i.e. c("GLM", "MLM", "FarmCPU")
# maxLine: when the priority is 'memory', users can change this parameter to limit the memory
# memo: a marker added on output file name
# P: a start p value for each SNP
# method.sub, method.sub.final: method used in substitution process
# method.bin: EMMA or FaSTLMM
# bin.size: window size in genome
# bin.selection: a vector, how many windows selected
# Prior: four columns, SNP name, Chr, Pos, P
# maxLoop: maximum number of iterations
# threshold.output: output GWAS results only for SNPs with p value lower than the threshold.output
# iteration.output: whether to output results for FarmCPU iterations
# p.threshold: if all p values in the 1st iteration are bigger than p.threshold, FarmCPU stops
# QTN.threshold: Only SNPs have a more significant p value than QTN.threshold have chance to be selected as pseudo QTNs
# bound: maximum number of SNPs selected as pseudo QTNs for each iteration
# outward: the direction of circular Manhattan plot
# permutation.threshold: if use a permutation cutoff or not (bonferroni cutoff)
# permutation.rep: number of permutation replicates
# col: for color of points in each chromosome on manhattan plot
# plot.type: "b" (both Manhattan plot and qq plot will be draw) or "q" (qq plot only)
# file.output: whether to output files or not
# file: figure formats, "jpg", "tiff"
# dpi: resolution
# threshold: a cutoff line on manhattan plot, 0.05/marker size
# Ncluster: number of colors used for drawing PC 1 and PC 2
# signal.cex: point size on output figures
#
# Output:
# results: a m * 2 matrix, the first column is the SNP effect, the second column is the P values
#
# Authors: Lilin Yin and Xiaolei Liu
# Build date: Aug 30, 2017
# Last update: January 27, 2017
###################################################################################################
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    linux <- R.ver == 'Linux'
    mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
    
    if(wind)	ncpus <- 1
    if(r.open && ncpus>1 && mac){
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    #if(r.open && ncpus>1 && !mac){
    #setMKLthreads(1)
    #}
    
    MVP.Version(TRUE)
    if(nrow(phe) != ncol(geno))	stop("The number of individuals in phenotype and genotype doesn't match!")
    #list -> matrix
    map <- as.matrix(map)
    na.index <- NULL
    if(!is.null(CV.GLM)){
        CV.GLM <- as.matrix(CV.GLM)
    if(nrow(CV.GLM) != ncol(geno))	stop("The number of individuals in covariates and genotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.GLM), arr.ind=T)[, 1])
    }
    if(!is.null(CV.MLM)){
        CV.MLM <- as.matrix(CV.MLM)
        if(nrow(CV.MLM) != ncol(geno))	stop("The number of individuals in covariates and genotype doesn't match!")
       na.index <- c(na.index, which(is.na(CV.MLM), arr.ind=T)[, 1])
    }
    if(!is.null(CV.FarmCPU)){
        CV.FarmCPU <- as.matrix(CV.FarmCPU)
        if(nrow(CV.FarmCPU) != ncol(geno))	stop("The number of individuals in covariates and genotype doesn't match!")
        na.index <- c(na.index, which(is.na(CV.FarmCPU), arr.ind=T)[, 1])
    }
    na.index <- unique(na.index)
    
    #remove samples with missing phenotype
    seqTaxa = which(!is.na(phe[,2]))
    if(length(na.index) != 0)	seqTaxa <- intersect(seqTaxa, c(1:nrow(phe))[-na.index])
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
    
    if(!is.null(nPC)){
        ipca <- MVP.PCA(M=geno, perc=perc, pcs.keep=nPC)$PCs
        if(file.output){
        filebck <- paste("MVP.", colnames(phe)[2], memo, ".pc.bin", sep="")
        filedes <- paste("MVP.", colnames(phe)[2], memo, ".pc.desc", sep="")
        if(file.exists(filebck) & file.exists(filedes)){
            PC.backed <- attach.big.matrix(filedes)
        }else{
                    PC.backed<-big.matrix(nrow(ipca), ncol(ipca), type="double", backingfile=filebck,
                descriptorfile=filedes)
        }
            PC.backed[, ] <- ipca[, ]
            flush(PC.backed)
            #print("Preparation for PC matrix is done!")
            rm(list=c("PC.backed"))
        }
        #CV for GLM
        if(glm.run){
            if(!is.null(CV.GLM)){
                if(!is.null(nPC.GLM)){
                    CV.GLM <- cbind(ipca[,1:nPC.GLM], CV.GLM)
                }
            }else{
                CV.GLM <- ipca[,1:nPC.GLM]
            }
        }
        
        #CV for MLM
        if(mlm.run){
            if(!is.null(CV.MLM)){
                if(!is.null(nPC.MLM)){
                    CV.MLM <- cbind(ipca[,1:nPC.MLM], CV.MLM)
                }
            }else{
                CV.MLM <- ipca[,1:nPC.MLM]
            }
        }
        
        #CV for FarmCPU
        if(farmcpu.run){
            if(!is.null(CV.FarmCPU)){
                if(!is.null(nPC.FarmCPU)){
                    CV.FarmCPU <- cbind(ipca[,1:nPC.FarmCPU], CV.FarmCPU)
                }
            }else{
                CV.FarmCPU <- ipca[,1:nPC.FarmCPU]
            }
        }
        
    }
  
    #GWAS
    print("GWAS Start...")
    
    if(glm.run){
        print("General Linear Model (GLM) Start ...")
        glm.results <- MVP.GLM(phe=phe, geno=geno, priority=priority, CV=CV.GLM, cpu=ncpus, memo="MVP.GLM", bar=bar);gc()
        colnames(glm.results) <- c("effect", paste(colnames(phe)[2],"GLM",sep="."))
        if(file.output)	write.csv(cbind(map,glm.results), paste("MVP.",colnames(phe)[2],".GLM", ".csv", sep=""), row.names=FALSE)
    }

    if(mlm.run){
        print("Mixed Linear Model (MLM) Start ...")
        mlm.results <- MVP.MLM(phe=phe, geno=geno, K=K, priority=priority, CV=CV.MLM, cpu=ncpus, bar=bar, maxLine=maxLine, vc.method=vc.method, file.output=file.output, memo="MVP.MLM");gc()
        colnames(mlm.results) <- c("effect", paste(colnames(phe)[2],"MLM",sep="."))
        if(file.output)	write.csv(cbind(map,mlm.results), paste("MVP.",colnames(phe)[2],".MLM", ".csv", sep=""), row.names=FALSE)
    }
    
    if(farmcpu.run){
        print("FarmCPU Start ...")
        farmcpu.results <- MVP.FarmCPU(phe=phe, geno=geno, map=map, priority=priority, CV=CV.FarmCPU, ncpus=ncpus, memo="MVP.FarmCPU", P=P, method.sub=method.sub, method.sub.final=method.sub.final, method.bin=method.bin, bin.size=bin.size, bin.selection=bin.selection, Prior=Prior, maxLoop=maxLoop, threshold.output=threshold.output, iteration.output=iteration.output, p.threshold=p.threshold, QTN.threshold=QTN.threshold, bound=NULL)
        colnames(farmcpu.results) <- c("effect", paste(colnames(phe)[2],"FarmCPU",sep="."))
        if(file.output)	write.csv(cbind(map,farmcpu.results), paste("MVP.",colnames(phe)[2],".FarmCPU", ".csv", sep=""), row.names=FALSE)
    }
    
    MVP.return <- list(map=map, glm.results=glm.results, mlm.results=mlm.results, farmcpu.results=farmcpu.results)
    
    if(permutation.threshold){
        set.seed(12345)
        i=1
            for(i in 1:permutation.rep){
            index = 1:nrow(phe)
            index.shuffle = sample(index,length(index),replace=F)
            myY.shuffle = phe
            myY.shuffle[,2] = myY.shuffle[index.shuffle,2]
            #GWAS using t.test...
            myPermutation = MVP.GLM(phe=myY.shuffle[,c(1,2)], geno=geno, cpu=ncpus)
            pvalue = min(myPermutation[,2],na.rm=T)
            if(i==1){
                    pvalue.final=pvalue
               }else{
                    pvalue.final=c(pvalue.final,pvalue)
            }
        }#end of permutation.rep
        permutation.cutoff = sort(pvalue.final)[ceiling(permutation.rep*0.05)]
        threshold = permutation.cutoff * m
    }
    print(paste("Significance Level: ", threshold/m, sep=""))
    if(file.output){
        print("Visualization Start...")
        print("Phenotype distribution Plotting...")
        MVP.Hist(phe=phe, file=file, col=col, dpi=dpi)
        #plot3D <- class(try(library("rgl"),silent=TRUE)) != "try-error"
        plot3D <- TRUE
        if(!is.null(nPC))	MVP.PCAplot(ipca[,1:3], col=col, plot3D=plot3D, Ncluster=Ncluster, file=file, dpi=dpi, box=box)
        MVP.Report(MVP.return, col=col, box=box, plot.type=c("c","m","q","d"), file.output=TRUE, file=file, dpi=dpi, threshold=threshold/m, signal.cex=signal.cex, outward=outward)
        if(sum(c(is.null(glm.results), is.null(mlm.results), is.null(farmcpu.results))) < 2)	MVP.Report(MVP.return, col=col, plot.type=c("m","q"), multracks=TRUE, outward=outward, file.output=TRUE, file=file, dpi=dpi, box=box, threshold=threshold/m, signal.cex=signal.cex)
    }
    MVP.Version(FALSE)
    return(MVP.return)
}#end of MVP function
