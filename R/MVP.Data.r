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
#' @param fileMVP Genotype in MVP format
#' @param fileVCF Genotype in VCF format
#' @param fileHMP Genotype in hapmap format
#' @param fileBed Genotype in PLINK binary format
#' @param fileInd Individual name file
#' @param fileNum Genotype in numeric format; pure 0, 1, 2 matrix; m * n, m is marker size, n is sample size
#' @param filePhe Phenotype, the first column is taxa name, the subsequent columns are traits
#' @param fileMap SNP map information, there are three columns, including SNP_ID, Chromosome, and Position
#' @param fileKin Kinship that represents relationship among individuals, n * n matrix, n is sample size
#' @param filePC Principal components, n*npc, n is sample size, npc is number of top columns of principal components
#' @param out prefix of output file name
#' @param sep.num seperator for numeric file.
#' @param auto_transpose Whether to automatically transpose numeric genotypes,
#' the default is True, which will identify the most one of the rows or 
#' columns as a marker, If set to False, the row represents the marker and
#' the column represents the individual.
#' @param sep.map seperator for map file.
#' @param sep.phe seperator for phenotype file.
#' @param sep.kin seperator for Kinship file.
#' @param sep.pc seperator for PC file.
#' @param type.geno type parameter in bigmemory, genotype data. The default is char, it is highly recommended *NOT* to modify this parameter.
#' @param pheno_cols Extract which columns of the phenotype file (including individual IDs)
#' @param SNP.impute "Left", "Middle", "Right", or NULL for skip impute.
#' @param maxLine number of SNPs, only used for saving memory when calculate kinship matrix
#' @param priority "speed" or "memory"
#' @param pcs.keep how many PCs to keep
#' @param verbose whether to print detail.
#' @param ncpus The number of threads used, if NULL, (logical core number - 1) is automatically used
#' @param ... Compatible with DEPRECATED parameters.
#'
#' @export
#' 
#' @return NULL
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' k.desc, k.bin: Kinship matrix in bigmemory format
#' pc.desc, pc.bin: PC matrix in bigmemory format
#' Requirement: fileHMP, fileBed, and fileNum can not input at the same time
#' @examples 
#' \donttest{
#' bfilePath <- file.path(system.file("extdata", "02_bfile", package = "rMVP"), "mvp")
#' opts <- options(rMVP.OutputLog2File = FALSE)
#' 
#' MVP.Data(fileBed=bfilePath, out=tempfile("outfile"), ncpus=1)
#' 
#' options(opts)
#' }
MVP.Data <- function(fileMVP = NULL, fileVCF = NULL, fileHMP = NULL, fileBed = NULL, fileNum = NULL, fileMap = NULL,
                     filePhe = NULL, fileInd = NULL, fileKin = NULL, filePC = NULL, out = "mvp", sep.num = "\t",
                     auto_transpose = TRUE, sep.map = "\t", sep.phe = "\t", sep.kin = "\t", sep.pc = "\t",
                     type.geno = "char", pheno_cols = NULL, SNP.impute = "Major", maxLine = 10000, priority = "speed",
                     pcs.keep = 5, verbose = TRUE, ncpus = NULL, ...) {
    
    logging.initialize("MVP.Data", dirname(out))
    
    logging.log("Preparing data for MVP...\n", verbose = verbose)
    
    # Parameter compatible upgrade
    params <- list(...)
    deprecated <- c(
        "sep.vcf",
        "vcf.jump",
        "SNP.effect",
        "sep.hmp",
        "type.kin",
        "type.pc",
        "type.map",
        "maxRecord",
        "maxInd"
    )
    for (x in deprecated) {
        if (x %in% names(params)) {
            warning(paste(x, "has been DEPRECATED."))
        }
    }
    
    # Check Data Input
    if (is.null(ncpus)) { ncpus <- detectCores() - 1 }
    
    geno_files <- !sapply(list(
        fileMVP, fileVCF, fileHMP, fileBed, fileNum, fileMap
        ), is.null)
    
    flag <- paste(sapply(strsplit(as.character(geno_files), ''), `[[`, 1), collapse = '')   # flag = 'TFFFFF'
    
    # convert genotype file
    error_input <- function(geno_files) {
        if (length(which(geno_files[1:5])) != 1) {
            stop("Please input only one genotype data format!")
        }
        
        if (length(which(geno_files[5:6])) == 1) {
            stop("Both Map and Numeric genotype files are needed!")
        }
    }
    switch(flag,
           # fileMVP, fileVCF, fileHMP, fileBed, fileNum, fileMap
           # TODO: Fix it
           # TFFFFF = MVP.Data.MVP2MVP(),
           FTFFFF = 
               MVP.Data.VCF2MVP(
                   vcf_file = fileVCF, 
                   out = out,
                   verbose = verbose,
                   threads = ncpus
               ),
           FFTFFF = 
               MVP.Data.Hapmap2MVP(
                   hmp_file = fileHMP, 
                   out = out,
                   verbose = verbose
               ),
           FFFTFF = 
               MVP.Data.Bfile2MVP(
                   bfile = fileBed, 
                   out = out, 
                   maxLine = maxLine, 
                   priority = priority, 
                   type.geno = type.geno,
                   verbose = verbose,
                   threads = ncpus
               ),
           FFFFTT = 
               MVP.Data.Numeric2MVP( 
                   num_file = fileNum,
                   map_file = fileMap,
                   out = out, 
                   maxLine = maxLine, 
                   priority = priority, 
                   type.geno = type.geno,
                   auto_transpose = auto_transpose,
                   verbose = verbose
               ),
           error_input(geno_files)
    )
    
    # phenotype
    if (!is.null(filePhe)) {
        MVP.Data.Pheno(
            pheno_file = filePhe, 
            out = out, 
            header = TRUE,
            cols = pheno_cols, 
            sep = sep.phe
            # , missing = missing
        )
    }
    
    # impute
    if (!is.null(SNP.impute)) {
        MVP.Data.impute(
            mvp_prefix = out,
            method = SNP.impute,
            ncpus = ncpus,
            verbose = verbose
        )
    }
    
    # get kin
    if (!is.null(fileKin) && fileKin != FALSE) {
        K <- MVP.Data.Kin(
            fileKin = fileKin, 
            mvp_prefix = out, 
            priority = priority, 
            sep = sep.kin,
            cpus = ncpus
        )
    }
    
    # get pc
    if (!is.null(filePC) && filePC != FALSE) {
        MVP.Data.PC(
            filePC = filePC, 
            mvp_prefix = out, 
            K = K[,],
            pcs.keep = pcs.keep,
            priority = priority, 
            sep = sep.pc,
            cpus = ncpus
        )
    }


    logging.log("MVP data prepration accomplished successfully!\n", verbose = verbose)
} # end of MVP.Data function

#' MVP.Data.VCF2MVP: To transform vcf data to MVP package
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' Accept the | or / separated markers, any variant sites that are not 0 or 1 will be considered NA.
#' 
#' @param vcf_file Genotype in VCF format
#' @param out the name of output file
#' @param maxLine the max number of line to write to big matrix for each loop
#' @param type.geno the type of genotype elements
#' @param threads number of thread for transforming
#' @param verbose whether to print the reminder
#'
#' @return number of individuals and markers.
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' @examples 
#' \donttest{
#' vcfPath <- system.file("extdata", "01_vcf", "mvp.vcf", package = "rMVP")
#' 
#' MVP.Data.VCF2MVP(vcfPath, tempfile("outfile"), threads=1)
#' }
#' 
MVP.Data.VCF2MVP <- function(vcf_file, out='mvp', maxLine = 1e4, type.geno='char', threads=1, verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    vcf_file <- normalizePath(vcf_file, mustWork = TRUE)
    # check old file
    backingfile <- paste0(basename(out), ".geno.bin")
    descriptorfile <- paste0(basename(out), ".geno.desc")
    remove_bigmatrix(out)
    
    # parser map
    logging.log("Reading file...\n", verbose = verbose)
    scan <- vcf_parser_map(vcf_file = vcf_file, out = out)
    m <- scan$m
    n <- scan$n
    logging.log(paste0("inds: ", n, "\tmarkers:", m, '\n'), verbose = verbose)
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = m,
        ncol = n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    vcf_parser_genotype(vcf_file = vcf_file, pBigMat = bigmat@address, maxLine = maxLine, threads = threads, verbose = verbose)
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
    return(invisible(c(m, n)))
}

#' MVP.Data.Bfile2MVP: To transform plink binary data to MVP package
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param bfile Genotype in binary format (.bed, .bim, .fam)
#' @param out the name of output file
#' @param maxLine the max number of line to write to big matrix for each loop
#' @param priority 'memory' or 'speed'
#' @param type.geno the type of genotype elements
#' @param threads number of thread for transforming
#' @param verbose whether to print the reminder
#'
#' @return number of individuals and markers.
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' @examples 
#' \donttest{
#' bfilePath <- file.path(system.file("extdata", "02_bfile", package = "rMVP"), "mvp")
#' 
#' MVP.Data.Bfile2MVP(bfilePath, tempfile("outfile"), threads=1)
#' }
#' 
MVP.Data.Bfile2MVP <- function(bfile, out='mvp', maxLine=1e4, priority='speed', type.geno='char', threads=0, verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    bim_file <- normalizePath(paste0(bfile, '.bim'), mustWork = TRUE)
    fam_file <- normalizePath(paste0(bfile, '.fam'), mustWork = TRUE)
    bed_file <- normalizePath(paste0(bfile, '.bed'), mustWork = TRUE)
    # check old file
    backingfile <- paste0(basename(out), ".geno.bin")
    descriptorfile <- paste0(basename(out), ".geno.desc")
    remove_bigmatrix(out)
    
    # parser map
    logging.log("Reading file...\n", verbose = verbose)
    m <- MVP.Data.Map(bim_file, out = out, cols = c(2, 1, 4, 6, 5), header = FALSE)
    
    # parser phenotype, ind file
    fam <- read.table(fam_file, header = FALSE)
    n <- nrow(fam)
    write.table(fam[, 2], paste0( out, '.geno.ind'), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    logging.log(paste0("inds: ", n, "\tmarkers:", m, '\n'), verbose = verbose)
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = m,
        ncol = n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    if (priority == "speed") { maxLine <- -1 }
    read_bfile(bed_file = bed_file, pBigMat = bigmat@address, maxLine = maxLine, threads = threads, verbose = verbose)
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
    return(invisible(c(m, n)))
}

#' MVP.Data.Hapmap2MVP: To transform Hapmap data to MVP package
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param hmp_file Genotype in Hapmap format
#' @param out the name of output file
#' @param maxLine the max number of line to write to big matrix for each loop
#' @param type.geno the type of genotype elements
#' @param threads number of thread for transforming
#' @param verbose whether to print the reminder
#'
#' @return number of individuals and markers.
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' @examples 
#' \donttest{
#' hapmapPath <- system.file("extdata", "03_hapmap", "mvp.hmp.txt", package = "rMVP")
#' 
#' MVP.Data.Hapmap2MVP(hapmapPath, tempfile("outfile"), threads=1)
#' }
#' 
MVP.Data.Hapmap2MVP <- function(hmp_file, out='mvp', maxLine = 1e4, type.geno='char', threads=1, verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    hmp_file <- normalizePath(hmp_file, mustWork = TRUE)
    # check old file
    backingfile <- paste0(basename(out), ".geno.bin")
    descriptorfile <- paste0(basename(out), ".geno.desc")
    remove_bigmatrix(out)
    
    # parser map
    logging.log("Reading file...\n", verbose = verbose)
    logging.log("Writing map into file\n", verbose = verbose)
    scan <- hapmap_parser_map(hmp_file[1], out)
    m <- scan$m
    n <- scan$n
    logging.log(paste0("inds: ", n, "\tmarkers:", m, '\n'), verbose = verbose)
    
    # parse genotype
    bigmat <- filebacked.big.matrix(
        nrow = m,
        ncol = n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    hapmap_parser_genotype(hmp_file = hmp_file, Major = scan$Major, pBigMat = bigmat@address, maxLine = maxLine, threads = threads, verbose = verbose)
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
    return(invisible(c(m, n)))
}

#' MVP.Data.Numeric2MVP: To transform Numeric data to MVP package
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param num_file Genotype in Numeric format (0,1,2)
#' @param map_file Genotype map file, SNP_name, Chr, Pos 
#' @param out the name of output file
#' @param maxLine the max number of line to write to big matrix for each loop
#' @param priority 'memory' or 'speed'
#' @param row_names whether the numeric genotype has row names
#' @param col_names whether the numeric genotype has column names
#' @param type.geno the type of genotype elements
#' @param auto_transpose whether to detecte the row and column
#' @param verbose whether to print the reminder
#'
#' @return number of individuals and markers.
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' @examples 
#' \donttest{
#' numericPath <- system.file("extdata", "04_numeric", "mvp.num", package = "rMVP")
#' mapPath <- system.file("extdata", "04_numeric", "mvp.map", package = "rMVP")
#' MVP.Data.Numeric2MVP(numericPath, mapPath, tempfile("outfile"))
#' }
#' 
MVP.Data.Numeric2MVP <- function(num_file, map_file, out='mvp', maxLine=1e4, priority='speed', row_names=FALSE, col_names=FALSE, type.geno='char', auto_transpose=TRUE, verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    num_file <- normalizePath(num_file, mustWork = TRUE)
    map_file <- normalizePath(map_file, mustWork = TRUE)
    # check old file
    backingfile <- paste0(basename(out), ".geno.bin")
    descriptorfile <- paste0(basename(out), ".geno.desc")
    remove_bigmatrix(out)
    
    # detecte n(ind) and m(marker)
    logging.log("Reading file...\n", verbose = verbose)
    scan <- numeric_scan(num_file)
    n <- scan$n
    m <- scan$m

    transposed <- FALSE
    if (auto_transpose & (m < n)) {
        message("WARNING: nrow < ncol detected, has been automatically transposed.")
        transposed <- TRUE
        t <- n; n <- m; m <- t;
    }
    logging.log(paste0("inds: ", n, "\tmarkers:", m, '\n'), verbose = verbose)
    
    # define bigmat
    bigmat <- filebacked.big.matrix(
        nrow = m,
        ncol = n,
        type = type.geno,
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    # convert to bigmat - speed
    if (priority == "speed") {
        opts <- options(bigmemory.typecast.warning = FALSE)
        on.exit(options(opts))
        
        # detecte sep
        con <- file(num_file, open = 'r')
        line <- readLines(con, 1)
        close(con)
        sep <- substr(line, 2, 2)
        
        # load geno
        suppressWarnings(
            geno <- read.big.matrix(num_file, header = FALSE, sep = sep)
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
        if (col_names) { readLines(con, n = 1) }
        while (TRUE) {
            line = readLines(con, n = maxLine)

            len <- length(line)
            if (len == 0) { break }

            line <- do.call(rbind, strsplit(line, '\\s+'))
            if (row_names) { line <- line[, 2:ncol(line)]}
            if (transposed) {
                bigmat[, (i + 1):(i + ncol(line))] <- line
                i <- i + ncol(line)
                percent <- 100 * i / n
            } else {
                bigmat[(i + 1):(i + nrow(line)), ] <- line
                i <- i + nrow(line)
                percent <- 100 * i / m
            }

            logging.log(paste0("Written into MVP File: ", percent, "%"), verbose = verbose)
        }
        close(con)
    }
    
    file.copy(map_file, paste0(out, ".geno.map"))
    
    flush(bigmat)
    gc()
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
    return(invisible(c(m, n)))
}

#' MVP.Data.MVP2Bfile: To transform MVP data to binary format
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param bigmat Genotype in bigmatrix format (0,1,2)
#' @param map the map file
#' @param pheno the phenotype file
#' @param out the name of output file
#' @param threads number of thread for transforming
#' @param verbose whether to print the reminder
#'
#' @return NULL
#' Output files:
#' .bed, .bim, .fam
#' 
#' @examples 
#' \donttest{
#' bigmat <- as.big.matrix(matrix(1:6, 3, 2))
#' map <- matrix(c("rs1", "rs2", "rs3", 1, 1, 1, 10, 20, 30), 3, 3)
#' 
#' MVP.Data.MVP2Bfile(bigmat, map, out=tempfile("outfile"), threads=1)
#' }
#' 
MVP.Data.MVP2Bfile <- function(bigmat, map, pheno=NULL, out='mvp.plink', threads=1, verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    # write bed file
    write_bfile(bigmat@address, out, threads=threads, verbose = verbose)
    
    # write fam
    #  1. Family ID ('FID')
    #  2. Within-family ID ('IID'; cannot be '0')
    #  3. Within-family ID of father ('0' if father isn't in dataset)
    #  4. Within-family ID of mother ('0' if mother isn't in dataset)
    #  5. Sex code ('1' = male, '2' = female, '0' = unknown)
    #  6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
    if (is.null(pheno)) {
        ind <- paste0("ind", 1:ncol(bigmat))
        pheno <- rep(-9, ncol(bigmat))
        message("pheno is NULL, automatically named individuals.")
    } else if (ncol(pheno) == 1) {
        ind <- pheno[, 1]
        pheno <- rep(-9, ncol(bigmat))
    } else {
        if (ncol(pheno) > 2) { 
            message("Only the first phenotype is written to the fam file, and the remaining ", ncol(pheno) - 1, " phenotypes are ignored.")
        }
        ind <- pheno[, 1]
        pheno <- pheno[, 2]
    }
    
    fam <- cbind(ind, ind, 0, 0, 0, pheno)
    write.table(fam, paste0(out, '.fam'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
    
    # write bim
    #  1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
    #  2. Variant identifier
    #  3. Position in morgans or centimorgans (safe to use dummy value of '0')
    #  4. Base-pair coordinate (normally 1-based, but 0 ok; limited to 231-2)
    #  5. Allele 1 (corresponding to clear bits in .bed; usually minor)
    #  6. Allele 2 (corresponding to set bits in .bed; usually major)
    bim <- cbind(map[, 2], map[, 1], 0, map[, 3], 0, 0)
    write.table(bim, paste0(out, '.bim'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
    t2 <- as.numeric(Sys.time())
    logging.log("Done within", format_time(t2 - t1), "\n", verbose = verbose)
}

#' MVP.Data.Pheno: To clean up phenotype file
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param pheno_file the name of phenotype file
#' @param out the name of output file
#' @param cols selected columns
#' @param header whether the file contains header
#' @param sep seperator of the file
#' @param missing the missing value
#' @param verbose whether to print detail.
#'
#' @return NULL
#' Output files:
#' cleaned phenotype file
#' @examples 
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' 
#' MVP.Data.Pheno(phePath, out=tempfile("outfile"))
#' 
MVP.Data.Pheno <- function(pheno_file, out='mvp', cols=NULL, header=TRUE, sep='\t', missing=c(NA, 'NA', '-9', 9999), verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    # read data
    if (!is.vector(pheno_file)) { pheno_file <- c(pheno_file) }

    # phenotype files
    phe <- read.delim(pheno_file, sep = sep, header = header)
    
    # auto select columns
    if (is.null(cols)) {
        cols <- c(1:ncol(phe))
    }
    
    # check phenotype file
    if (length(cols) < 2) {
        stop("ERROR: At least 2 columns in the phenotype file should be specified, please check the parameter 'sep.phe'. ")
    }
    phe[, cols[1]] <- sapply(phe[, cols[1]], function(x){gsub("^\\s+|\\s+$", "", x)}) 
    
    # read geno ind list
    geno.ind.file <- paste0(out, '.geno.ind')
    if (file.exists(geno.ind.file)) {
        # read from file
        geno.ind <- read.table(geno.ind.file, stringsAsFactors = FALSE)[, 1]
        overlap.ind <- intersect(geno.ind, phe[, cols[1]])
        if (length(overlap.ind) == 0) {
            logging.log(paste0("Phenotype individuals: ", paste(phe[, cols[1]][1:5], collapse = ", "), "..."), "\n", verbose = verbose)
            logging.log(paste0("Genotype individuals: ", paste(geno.ind[1:5], collapse = ", "), "..."), "\n", verbose = verbose)
            stop("No common individuals between phenotype and genotype!")
        } else {
            logging.log(paste(length(overlap.ind), "common individuals between phenotype and genotype."), "\n", verbose = verbose)
        }
    } else {
        # use ind. name from phenotypefile
        geno.ind <- phe[, cols[1]]
    }
    
    # merge
    pheno <- phe[match(geno.ind, phe[, cols[1]]), cols]
    pheno[1] <- geno.ind
    
    # rename header
    colnames(pheno)[1] <- 'Taxa'
    if (!header)  {
        traits <- 2:ncol(pheno)
        colnames(pheno)[traits] <- paste0('t', traits - 1)
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
    write.table(pheno, paste0(out, '.phe'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for PHENOTYPE data is Done within", format_time(t2 - t1), "\n", verbose = verbose)
}

#' MVP.Data.Map: To check map file
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param map the name of map file or map object(data.frame or matrix)
#' @param out the name of output file
#' @param cols selected columns
#' @param header whether the file contains header
#' @param sep seperator of the file
#' @param verbose whether to print detail.
#' 
#' @return 
#' Output file:
#' <out>.map
#' 
#' @examples 
#' mapPath <- system.file("extdata", "05_mvp", "mvp.geno.map", package = "rMVP")
#' 
#' MVP.Data.Map(mapPath, tempfile("outfile"))
#' 
MVP.Data.Map <- function(map, out='mvp', cols=1:5, header=TRUE, sep='\t', verbose=TRUE) {
    t1 <- as.numeric(Sys.time())
    if (is.character(map) && !is.data.frame(map)) {
        map <- read.table(map, header = header, stringsAsFactors = FALSE)
    }
    map <- map[, cols]
    colnames(map) <- c("SNP", "CHROM", "POS", "REF", "ALT")
    if (length(unique(map[, 1])) != nrow(map)) {
        warning("WARNING: SNP is not unique and has been automatically renamed.")
        map[, 1] <- paste(map[, 2], map[, 3], sep = "-")
    }
    allels <- map[, 4:5]
    allels[allels == 0] <- '.'
    map[, 4:5] <- allels
    
    write.table(map, paste0(out, ".geno.map"), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    t2 <- as.numeric(Sys.time())
    logging.log("Preparation for MAP data is done within", format_time(t2 - t1), "\n", verbose = verbose)
    return(nrow(map))
}


#' Principal component analysis
#'
#' @param filePC Principal components, n*npc, n is sample size, npc is number of top columns of principal components
#' @param mvp_prefix Prefix for mvp format files
#' @param K Kinship matrix
#' @param out prefix of output file name
#' @param pcs.keep how many PCs to keep
#' @param priority speed or memory
#' @param sep seperator for PC file.
#' @param cpus the number of cpu
#' @param verbose whether to print detail.
#' 
#' @export
#' @return 
#' Output file:
#' <out>.pc.bin
#' <out>.pc.desc
#' @examples
#' \donttest{
#' geno <- file.path(system.file("extdata", "06_mvp-impute", package = "rMVP"), "mvp.imp")
#' 
#' MVP.Data.PC(TRUE, mvp_prefix=geno, out=tempfile("outfile"), cpus=1)
#' }
#' 
MVP.Data.PC <- function(
    filePC=TRUE, 
    mvp_prefix='mvp', 
    K=NULL, 
    out=NULL,  
    pcs.keep=5,
    priority='speed',
    sep='\t',
    cpus=1,
    verbose=TRUE
){
    if (is.null(out)) out <- mvp_prefix
    
    # check old file
    backingfile <- paste0(basename(out), ".pc.bin")
    descriptorfile <- paste0(basename(out), ".pc.desc")
    remove_bigmatrix(out, desc_suffix = ".pc.desc", bin_suffix = ".pc.bin")
    
    if (is.character(filePC)) {
        myPC <- read.big.matrix(filePC, header = FALSE, type = 'double', sep = sep)
    } else if (filePC == TRUE) {
        if(is.null(K)){
            geno <- attach.big.matrix(paste0(mvp_prefix, ".geno.desc"))
            myPC <- MVP.PCA(M=geno, pcs.keep = pcs.keep, priority=priority, cpu=cpus)
        }else{
            myPC <- MVP.PCA(K=K, pcs.keep = pcs.keep, priority=priority, cpu=cpus)
        }
    } else {
        stop("ERROR: The value of filePC is invalid.")
    }
    
    # define bigmat
    PC <- filebacked.big.matrix(
        nrow = nrow(myPC),
        ncol = ncol(myPC),
        type = 'double',
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    PC[, ] <- myPC[, ]
    flush(PC)
    logging.log("Preparation for PC matrix is done!", "\n", verbose = verbose)
}

#' Kinship
#'
#' @param fileKin Kinship that represents relationship among individuals, n * n matrix, n is sample size
#' @param mvp_prefix Prefix for mvp format files
#' @param out prefix of output file name
#' @param priority "speed" or "memory"
#' @param sep seperator for Kinship file.
#' @param cpus the number of cpu
#' @param verbose whether to print detail.
#'
#' @export
#' @return 
#' Output file:
#' <out>.kin.bin
#' <out>.kin.desc
#' @examples
#' \donttest{
#' geno <- file.path(system.file("extdata", "06_mvp-impute", package = "rMVP"), "mvp.imp")
#' 
#' MVP.Data.Kin(TRUE, mvp_prefix=geno, out=tempfile("outfile"), cpus=1)
#' }
#' 
MVP.Data.Kin <- function(
    fileKin=TRUE, 
    mvp_prefix='mvp', 
    out=NULL, 
    priority='speed', 
    sep='\t',
    cpus=1,
    verbose=TRUE
) {
    if (is.null(out)) out <- mvp_prefix
    
    # check old file
    backingfile <- paste0(basename(out), ".kin.bin")
    descriptorfile <- paste0(basename(out), ".kin.desc")
    remove_bigmatrix(out, desc_suffix = ".kin.desc", bin_suffix = ".kin.bin")
    
    if (is.character(fileKin)) {
        myKin <- read.big.matrix(fileKin, header = FALSE, type = 'double', sep = sep)
    } else if (fileKin == TRUE) {
        geno <- attach.big.matrix(paste0(mvp_prefix, ".geno.desc"))
        logging.log("Calculate KINSHIP using Vanraden method...", "\n", verbose = verbose)
        myKin <- MVP.K.VanRaden(geno, priority = priority, cpu = cpus)
    } else {
        stop("ERROR: The value of fileKin is invalid.")
    }
    
    # define bigmat
    Kinship <- filebacked.big.matrix(
        nrow = nrow(myKin),
        ncol = ncol(myKin),
        type = 'double',
        backingfile = backingfile,
        backingpath = dirname(out),
        descriptorfile = descriptorfile,
        dimnames = c(NULL, NULL)
    )
    
    Kinship[, ] <- myKin[, ]
    flush(Kinship)
    logging.log("Preparation for Kinship matrix is done!", "\n", verbose = verbose)
    return(Kinship)
}


#' MVP.Data.impute: To impute the missing genotype
#' Author: Haohao Zhang
#' Build date: Sep 12, 2018
#' 
#' @param mvp_prefix the prefix of mvp file
#' @param out the prefix of output file
#' @param method 'Major', 'Minor', "Middle"
#' @param ncpus number of threads for imputing
#' @param verbose whether to print the reminder
#' @return NULL
#' Output files:
#' imputed genotype file
#' @examples 
#' mvpPath <- file.path(system.file("extdata", "05_mvp", package = "rMVP"), "mvp")
#' 
#' MVP.Data.impute(mvpPath, tempfile("outfile"), ncpus=1)
#' 
MVP.Data.impute <- function(mvp_prefix, out=NULL, method='Major', ncpus=NULL, verbose=TRUE) {
    # input
    desc <- paste0(mvp_prefix, ".geno.desc")
    bigmat <- attach.big.matrix(desc)
    
    if (!hasNA(bigmat@address)) {
        message("No NA in genotype, imputation has been skipped.")
        return()
    }
    
    logging.log("Imputing...\n", verbose = verbose)
    
    opts <- options(bigmemory.typecast.warning = FALSE)
    on.exit(options(opts))
    
    if (is.null(ncpus)) ncpus <- detectCores()
    
    if (is.null(out)) {
        message("out is NULL, impute inplace.")
        outmat <- attach.big.matrix(desc)
    } else {
        # output to new genotype file.
        backingfile <- paste0(basename(out), ".geno.bin")
        descriptorfile <- paste0(basename(out), ".geno.desc")
        remove_bigmatrix(out)
        
        outmat <- filebacked.big.matrix(
            nrow = nrow(bigmat),
            ncol = ncol(bigmat),
            type = typeof(bigmat),
            backingfile = backingfile,
            backingpath = dirname(out),
            descriptorfile = descriptorfile,
            dimnames = c(NULL, NULL)
        )
        outmat[, ] <- bigmat[, ]
        file.copy(paste0(mvp_prefix, ".geno.ind"), paste0(out, ".geno.ind"))
        file.copy(paste0(mvp_prefix, ".geno.map"), paste0(out, ".geno.map"))
    }
    
    impute_marker(outmat@address, threads = ncpus, verbose = verbose)
    
    logging.log("Impute Genotype File is done!\n", verbose = verbose)
}

