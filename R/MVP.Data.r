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


#' MVP.Data: To prepare data for MVP package
#' Author: Xiaolei Liu, Lilin Yin and Haohao Zhang
#' Build date: Aug 30, 2016
#' Last update: Sep 12, 2018
#' 
#' @param fileVCF Genotype in VCF format
#' @param fileHMP Genotype in hapmap format
#' @param fileBed Genotype in PLINK binary format
#' @param fileNum Genotype in numeric format; pure 0, 1, 2 matrix; m * n, m is marker size, n is sample size
#' @param filePhe Phenotype, two columns, the first column is taxa name, the second column is trait
#' @param fileMap SNP map information, three columns: SNP name, Chr, Pos
#' @param fileKin Kinship, n * n matrix, n is sample size
#' @param filePC 
#' @param out A marker on output file name
#' @param sep.vcf seperator for hapmap, numeric, map, phenotype, kinship and PC files, respectively
#' @param sep.hmp 
#' @param sep.num 
#' @param sep.map 
#' @param sep.phe 
#' @param sep.kin 
#' @param sep.pc 
#' @param vcf.jump 
#' @param type.geno type parameter in bigmemory, genotype data
#' @param type.kin type parameter in bigmemory, Kinship
#' @param type.pc type parameter in bigmemory, PC
#' @param type.map type parameter in bigmemory, PC
#' @param SNP.effect "Add" or "Dom"
#' @param SNP.impute "Left", "Middle", "Right"
#' @param maxLine number of SNPs, only used for saving memory when calculate kinship matrix
#' @param maxRecord maximum number for markers
#' @param maxInd maximum number for individuals
#' @param priority 
#' @param perc 
#' @param pcs.keep 
#'
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' k.desc, k.bin: Kinship matrix in bigmemory format
#' pc.desc, pc.bin: PC matrix in bigmemory format
#' Requirement: fileHMP, fileBed, and fileNum can not input at the same time

MVP.Data <- function(fileMVP = NULL, fileVCF = NULL, fileHMP = NULL, fileBed = NULL, fileNum = NULL, fileMap = NULL,
                     filePhe = NULL, fileInd = NULL, fileKin = TRUE, filePC = TRUE, out = "mvp", sep.num = "\t", auto_transpose = T,
                     sep.map = "\t", sep.phe = "\t", sep.kin = "\t", sep.pc = "\t", type.geno = "char", pheno_cols = NULL,
                     SNP.impute = "Major", maxLine = 10000, priority = "speed", perc = 1, pcs.keep = 5, ...) {
    
    print("Preparing data for MVP...")
    
    # Parameter compatible upgrade
    params <- list(...)
    if ("sep.vcf" %in% names(params)) { message("WARNING: 'sep.vcf' has been DEPRECATED.") }
    if ("vcf.jump" %in% names(params)) { message("WARNING: 'vcf.jump' has been DEPRECATED.") }
    if ("SNP.effect" %in% names(params)) { message("WARNING: 'SNP.effect' has been DEPRECATED.") }
    if ("sep.hmp" %in% names(params)) { message("WARNING: 'sep.hmp' has been DEPRECATED.") }
    if ("type.kin" %in% names(params)) { message("WARNING: 'type.kin' has been DEPRECATED.") }
    if ("type.pc" %in% names(params)) { message("WARNING: 'type.pc' has been DEPRECATED.") }
    if ("type.map" %in% names(params)) { message("WARNING: 'type.map' has been DEPRECATED.") }
    if ("maxRecord" %in% names(params)) { message("WARNING: 'maxRecord' has been DEPRECATED.") }
    if ("maxInd" %in% names(params)) { message("WARNING: 'maxInd' has been DEPRECATED.") }
    
    # Check Data Input
    geno_files <- !sapply(list(
        fileMVP, fileVCF, fileHMP, fileBed, fileNum, fileMap
        ), is.null)
    
    flag <- paste(sapply(strsplit(as.character(geno_files), ''), `[[`, 1), collapse = '')
    
    error_input <- function(geno_files) {
        if (length(which(geno_files[1:5])) != 1) {
            stop("Please input only one genotype data format!")
        }
        
        if (length(which(geno_files[5:6])) == 1) {
            stop("Both Map and Numeric genotype files are needed!")
        }
    }
    
    # convert genotype file
    switch(flag,
           # fileMVP, fileVCF, fileHMP, fileBed, fileNum, fileMap
           TFFFFF = MVP.Data.MVP2MVP(),
           FTFFFF = MVP.Data.VCF2MVP(fileVCF, out),
           FFTFFF = MVP.Data.Hapmap2MVP(fileHMP, out),
           FFFTFF = MVP.Data.Bfile2MVP(fileBed, out, maxLine, priority, type.geno),
           FFFFTT = MVP.Data.Numeric2MVP(fileNum, out, maxLine, priority, type.geno, auto_transpose),
           error_input(geno_files)
    )
    print("Preparation for Genotype File is done!")
    
    # phenotype
    if (!is.null(filePhe)) {
        MVP.Data.Pheno(filePhe, out, pheno_cols, sep = sep.phe)
    }
    # impute
    MVP.Data.impute(joint(out, '.geno.desc'), SNP.impute)
    
    # get pc
    MVP.Data.PC(filePC, out, perc, pcs.keep, sep.pc)
    
    # get kin
    MVP.Data.Kin(fileKin, out, maxLine, priority, sep.kin)

    print("MVP data prepration accomplished successfully!")
}# end of MVP.Data function

joint <- function(...) { paste(..., sep = "") }

MVP.Data.VCF2MVP <- function(vcf_file, out, type.geno='char') {
    # check old file
    backingfile <- joint(out, ".geno.bin")
    descriptorfile <- joint(out, ".geno.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # parser map
    m_res <- vcf_parser_map(vcf_file, out)
    cat(joint("inds: ", m_res$n, "\tmarkers:", m_res$m, '\n'))
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = m_res$m,
        ncol = m_res$n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    vcf_parser_genotype(vcf_file, bigmat@address, m_res$pos)
}

MVP.Data.Bfile2MVP <- function(bfile, out='mvp', maxLine=1e4, priority='speed', type.geno='char') {
    # check old file
    backingfile <- joint(out, ".geno.bin")
    descriptorfile <- joint(out, ".geno.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # parser map
    map <- MVP.Data.Map(joint(bfile, '.bim'), out, c(2, 1, 4))
    
    # parser phenotype
    fam <- read.delim(joint(bfile, '.fam'), header = F)
    
    cat(joint("inds: ", nrow(fam), "\tmarkers:", nrow(map), '\n'))
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = nrow(map),
        ncol = nrow(fam),
        type = type.geno,
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    if (priority == "speed") {maxLine <- -1}
    read_bfile(bfile, bigmat@address, maxLine)
}

MVP.Data.Hapmap2MVP <- function(hapmap_file, out='mvp', type.geno='char') {
    # check old file
    backingfile <- joint(out, ".geno.bin")
    descriptorfile <- joint(out, ".geno.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # parser map
    m_res <- hapmap_parser_map(hapmap_file, out)
    cat(joint("inds: ", m_res$n, "\tmarkers:", m_res$m, '\n'))
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = m_res$m,
        ncol = m_res$n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    hapmap_parser_genotype(hapmap_file, bigmat@address)
}

MVP.Data.Numeric2MVP <- function(num_file, out='mvp', maxLine=1e4, priority='speed', type.geno='char', auto_transpose=T) {
    # check old file
    backingfile <- joint(out, ".geno.bin")
    descriptorfile <- joint(out, ".geno.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # detecte n(ind) and m(marker)
    scan <- numeric_scan(num_file)
    n <- scan$n
    m <- scan$m

    transposed <- FALSE
    if (auto_transpose & (m < n)) {
        message("WARNING: nrow < ncol detected, has been automatically transposed.")
        transposed <- TRUE
        t <- n; n <- m; m <- t;
    }
    cat(joint("inds: ", n, "\tmarkers:", m, '\n'))
    
    # define bigmat
    bigmat <- filebacked.big.matrix(
        nrow = m,
        ncol = n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    # convert to bigmat - speed
    if (priority == "speed") {
        options(bigmemory.typecast.warning = FALSE)
        
        # detecte sep
        con <- file(num_file, open = 'r')
        line <- readLines(con, 1)
        close(con)
        sep <- substr(line, 2, 2)
        
        # load geno
        suppressWarnings(
            geno <- read.big.matrix(num_file, head = FALSE, sep = sep) 
        )
        if (transposed) {
            bigmat[, ] <- t(geno[, ])
        } else {
            bigmat[, ] <- geno[, ]
        }
        rm("geno")
    }
    
    # convert to bigmat - memory
    if (priority == "memory") {
        i <- 0
        con <- file(num_file, open = 'r')
        if (col_name) { readLines(con, n = 1) }
        while (TRUE) {
            line = readLines(con, n = maxLine)

            len <- length(line)
            if (len == 0) { break }

            line <- do.call(rbind, strsplit(line, '\\s+'))
            if (row_name) { line <- line[, 2:ncol(line)]}
            if (transposed) {
                bigmat[, (i + 1):(i + length(line))] <- line
                i <- i + length(line)
                percent <- 100 * i / n
            } else {
                bigmat[(i + 1):(i + length(line)), ] <- line
                i <- i + length(line)
                percent <- 100 * i / m
            }

            cat(joint("Written into MVP File: ", percent, "%"))
        }
        close(con)
    }
    
    flush(bigmat)
    gc()
}

MVP.Data.Pheno <- function(pheno_file, out='mvp', cols=NULL, header=T, sep='\t', missing=c(NA, 'NA', '-9')) {
    # read data
    if (!is.vector(pheno_file)) { pheno_file <- c(pheno_file) }

    # phenotype files
    f <- read.delim(pheno_file, sep = sep, header = header)
    
    # auto select columns
    if (is.null(cols)) {
        c <- c(1:ncol(f))
    } else {
        c <- cols
    }
    
    if (length(c) < 2) {
        stop("ERROR: At least 2 columns in the phenotype file should be specified.")
    }
    
    # read geno ind list
    geno.ind.file <- joint(out, '.geno.ind')
    if (file.exists(geno.ind.file)) {
        geno.ind <- read.table(geno.ind.file)
    } else {
        geno.ind <- f[, c[1]]
    }
    
    # merge
    pheno <- merge(geno.ind, f[, c],  by = 1, all.x = T)

    
    # rename header
    colnames(pheno)[1] <- 'Taxa'
    if (!header)  {
        traits <- 2:ncol(pheno)
        colnames(pheno)[traits] <- joint('t', traits - 1)
    }
    
    # drop empty traits
    pheno[pheno %in% missing] <- NA
    drop = c()
    for (i in 2:ncol(pheno)) {
        if (all(is.na(pheno[, i]))) {
            drop = c(drop, i)
        }
    }
    if (length(drop) > 0) {
        pheno <- pheno[, -drop]
    }
    
    # Output
    write.table(pheno, joint(out, '.phe'), quote = F, sep = "\t", row.names = F, col.names = T)
    
    print("Preparation for PHENOTYPE data is done!")
    return(pheno)
}

MVP.Data.Map <- function(map_file, out='mvp', cols=c(1, 2, 3), header=T, sep='\t') {
    map <- read.table(map_file, header = header)
    map <- map[, cols]
    colnames(map) <- c("SNP", "CHROM", "POS")
    if (length(unique(map[, 1])) != nrow(map)) {
        warning("WARNING: SNP is not unique and has been automatically renamed.")
        map[, 1] <- apply(map[, c(2, 3)], 1, paste, collapse = "-")
    }
    write.table(map, joint(out, ".map"), row.names = F, col.names = T, sep = '\t', quote = F)
    return(map)
}

MVP.Data.PC <- function(filePC, mvp_prefix='mvp', out=NULL, perc=1, pcs.keep=5, sep='\t') {
    if (is.null(out)) out <- mvp_prefix
    
    # check old file
    backingfile <- joint(out, ".pc.bin")
    descriptorfile <- joint(out, ".pc.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # get pc
    if (is.character(filePC)) {
        myPC <- read.big.matrix(filePC, head = FALSE, type = 'double', sep = sep)
    } else if (filePC == TRUE) {
        geno <- attach.big.matrix(joint(mvp_prefix, ".geno.desc"))
        myPC <- MVP.PCA(geno, perc = perc, pcs.keep = pcs.keep)$PCs
    } else if (filePC == FALSE || is.null(filePC)) {
        return()
    } else {
        stop("ERROR: The value of filePC is invalid.")
    }
    
    # define bigmat
    PC <- filebacked.big.matrix(
        nrow = nrow(myPC),
        ncol = ncol(myPC),
        type = 'double',
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    PC[, ] <- myPC[, ]
    flush(PC)
    print("Preparation for PC matrix is done!")
}

MVP.Data.Kin <- function(fileKin, mvp_prefix='mvp', out=NULL, maxLine=1e4, priority='speed', sep='\t') {
    if (is.null(out)) out <- mvp_prefix
    
    # check old file
    backingfile <- joint(out, ".kin.bin")
    descriptorfile <- joint(out, ".kin.desc")
    if (file.exists(backingfile)) file.remove(backingfile)
    if (file.exists(descriptorfile)) file.remove(descriptorfile)
    
    # get kin
    if (is.character(fileKin)) {
        myKin <- read.big.matrix(fileKin, head = F, type = 'double', sep = sep)
    } else if (fileKin == TRUE) {
        geno <- attach.big.matrix(joint(mvp_prefix, ".geno.desc"))
        print("Calculate KINSHIP using Vanraden method...")
        myKin <- MVP.K.VanRaden(geno, priority = priority, maxLine = maxLine)
    } else if (fileKin == FALSE || is.null(fileKin)) {
        return()
    } else {
        stop("ERROR: The value of fileKin is invalid.")
    }
    
    # define bigmat
    Kinship <- filebacked.big.matrix(
        nrow = nrow(myKin),
        ncol = ncol(myKin),
        type = 'double',
        backingfile = backingfile,
        backingpath = ".",
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    Kinship[, ] <- myKin[, ]
    flush(Kinship)
    print("Preparation for Kinship matrix is done!")
}

# TODO: Very slow (inds: 6, markers:50703 ~ 30s @haohao's mbp)
MVP.Data.impute <- function(mvp_file, method = 'Major') {
    bigmat  <- attach.big.matrix(mvp_file)
    options(bigmemory.typecast.warning = FALSE)
    
    for (i in 1:nrow(bigmat)) {
        # get frequency 
        c <- table(bigmat[i, ])
        
        # get Minor / Major / Middle Gene
        if (method == 'Middle' | length(c) == 0) { A <- 1 }
        else if (method == 'Major') { A <- as.numeric(names(c[c == max(c)])) }
        else if (method == 'Minor') { A <- as.numeric(names(c[c == min(c)])) }
        
        # impute
        if (length(A) > 1) { A <- sample(A, 1) }
        bigmat[i, is.na(bigmat[i, ])] <- A
    }
    print("Impute Genotype File is done!")
    
    # biganalytics::apply(bigmat, 1, impute.marker, MISSING = MISSING, method = method)
}

MVP.Data.QC <- function() {
    
}

HMP.Numeralization <-
function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE){
    #Object: To convert character SNP genotpe to numerical
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang
    #Modified by Xiaolei Liu
    #Last update: May 20, 2017
    ##############################################################################################
    if(bit==1)  {
        x[x=="X"]="N"
        x[x=="-"]="N"
        x[x=="+"]="N"
        x[x=="/"]="N"
        x[x=="K"]="Z"
    }
    
    if(bit==2)  {
		#replace missing allels
		x[x=="XX"]="N"
		x[x=="--"]="N"
		x[x=="++"]="N"
		x[x=="//"]="N"
		x[x=="NN"]="N"
		x[x=="00"]="N"
		
		#replace false allels
		x[x=="CA"]="AC"
		x[x=="GA"]="AG"
		x[x=="TA"]="AT"
		x[x=="GC"]="CG"
		x[x=="TC"]="CT"
		x[x=="TG"]="GT"
    }
    
    n=length(x)
    lev=levels(as.factor(x))
    lev=setdiff(lev,"N")
    len=length(lev)

    #Genotype counts
    count=1:len
    for(i in 1:len){
        count[i]=length(x[(x==lev[i])])
    }
    
    if(Major.allele.zero){
        if(len>1 & len<=3){
            #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
            if(bit==1){
                count.temp = cbind(count, seq(1:len))
                if(len==3) count.temp = count.temp[-3,]
                count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
                if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
            }
            
            #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
            if(bit==2){
                count.temp = cbind(count, seq(1:len))
                if(len==3) count.temp = count.temp[-2,]
                count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
                if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
            }
            
            count = count[order]
            lev = lev[order]
            
        }   #End  if(len<=1 | len> 3)
    } #End  if(Major.allele.zero)
    
    #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
    if(bit==1 & len==3){
        temp=count[2]
        count[2]=count[3]
        count[3]=temp
    }
    position=order(count)
    
    #1status other than 2 or 3
    if(len<=1 | len> 3)x=0
    
    #2 status
    if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,2))
    
    #3 status
    if(bit==1){
        if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
    }else{
        if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],2,1)))
    }
    
    #missing data imputation
    if(impute=="Middle") {x[is.na(x)]=1 }
    
    if(len==3){
        if(impute=="Minor")  {x[is.na(x)]=position[1]-1}
        if(impute=="Major")  {x[is.na(x)]=position[len]-1}
    }else{
        if(impute=="Minor")  {x[is.na(x)]=2*(position[1]-1)}
        if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
    }
    
    #alternative genetic models
    if(effect=="Dom") x=ifelse(x==1,1,0)
    if(effect=="Left") x[x==1]=0
    if(effect=="Right") x[x==1]=2
    
    result=matrix(x,n,1)

    return(result)
}
