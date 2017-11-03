MVP.Data <-
function(fileHMP=NULL, fileBed=NULL, fileNum=NULL, filePhe=NULL, fileMap=NULL, fileKin=TRUE, filePC=TRUE, out=NULL, sep.hmp="\t", sep.num="\t", sep.map="\t", sep.phe="\t", sep.kin="\t", sep.pc="\t", type.geno="char", type.kin="double", type.pc="double", type.map="integer", SNP.effect="Add",
SNP.impute="Middle", maxLine=10000, maxRecord=1e9, maxInd=1e9, priority="speed", perc=1, pcs.keep=5){
##########################################################################################################
# Object: To prepare data for MVP package
#
# Input:
# fileHMP: Genotype in hapmap format
# fileBed: Genotype in PLINK binary format
# fileNum: Genotype in numeric format; pure 0, 1, 2 matrix; m * n, m is marker size, n is sample size
# filePhe: Phenotype, two columns, the first column is taxa name, the second column is trait
# fileMap: SNP map information, three columns: SNP name, Chr, Pos
# fileKin: Kinship, n * n matrix, n is sample size
# out: A marker on output file name
# sep.hmp, sep.num, sep.map, sep.phe, sep.kin, sep.pc: seperator for hapmap, numeric, map, phenotype, kinship and PC files, respectively
# type.geno: type parameter in bigmemory, genotype data
# type.kin: type parameter in bigmemory, Kinship
# type.pc: type parameter in bigmemory, PC
# type.map: type parameter in bigmemory, map
# SNP.effect: "Add" or "Dom"
# SNP.impute: "Left", "Middle", "Right"
# maxLine: number of SNPs, only used for saving memory when calculate kinship matrix
# maxRecord: maximum number for markers
# maxInd: maximum number for individuals
#
# Output files:
# genotype.desc, genotype.bin: genotype file in bigmemory format
# phenotype.phe: ordered phenotype file, same taxa order with genotype file
# map.map: SNP information
# k.desc, k.bin: Kinship matrix in bigmemory format
# pc.desc, pc.bin: PC matrix in bigmemory format
# Requirement: fileHMP, fileBed, and fileNum can not input at the same time
#
# Authors: Xiaolei Liu and Lilin Yin
# Build date: Aug 30, 2016
# Last update: May 25, 2017
##############################################################################################
    
    print("Preparing data for MVP...")
    
    #if(is.null(fileHMP)&is.null(fileNum))
    #stop("Hapmap or Numeric genotype is needed.")
    if(!is.null(fileHMP)&!is.null(fileNum)&!is.null(fileBed))
    stop("Please don't input more than one data format!")
    if((!is.null(fileNum) & is.null(fileMap)) | (is.null(fileNum) & !is.null(fileMap)))
    stop("Both Map and Numeric genotype files are needed!")
    if(is.null(out)) out="MVP"
    
    if(!is.null(filePhe)){
        myY <- read.delim(filePhe, sep=sep.phe, head = TRUE)
        taxa.y <- as.vector(myY[, 1])
        taxa <- taxa.y
    }
    
    if(!is.null(fileBed)){
        phe <- read.table(paste(fileBed, ".fam", sep=""), head=FALSE)[, -c(1,3:5)]
        write.table(phe, paste(out, ".phe", sep=""), row.names=FALSE, col.names=TRUE, sep=sep.phe, quote=FALSE)
        map <- read.table(paste(fileBed, ".bim", sep=""), head=FALSE)
        if(length(unique(map[, 2])) != nrow(map)){
			stop("The names of all SNPs must be unique, please check 'BIM' file!")
		}
        map <- map[, c(2, 1, 4)]
        colnames(map) <- c("SNP", "chrom", "pos")
        write.table(map, paste(out, ".map", sep=""), row.names=FALSE, col.names=TRUE, sep=sep.map, quote=FALSE)
        print("Reading binary file...")
        geno <- read.plink(paste(fileBed, ".bed", sep=""))[[1]]
        print("Reading binary files is done!")
        #geno <- t(geno); gc()
        #geno=as.data.frame(geno)
        bck <- paste(out, ".geno.bin", sep="")
        dsc <- paste(out, ".geno.desc", sep="")
        nmarkers <- ncol(geno)
        ns <- nrow(geno)
        myGeno.backed<-big.matrix(nmarkers, ns, type="char",
        backingfile=bck, descriptorfile=dsc)
        options(bigmemory.typecast.warning=FALSE)
        #print("Output BIG genotype...")
	inGENOFile=TRUE
		i <- 0
		printN <- unlist(strsplit(x=as.character(nmarkers), split="", fixed=TRUE))
		printIndex <- seq(0, (as.numeric(printN[1]) + 1) * 10^(length(printN)), 1000)[-1]
		Num.fun <- function(x){
			x <- data.matrix(as.data.frame(x))
			x[x==0]=2
			return(x)
		}
		while(inGENOFile){
			i <- i + maxLine
			if(i >= nmarkers){
				xx <- nmarkers
				inGENOFile <- FALSE
			}else{
				xx <- i
			}
			if(sum(i >= printIndex )!=0){
				printIndex <- printIndex[printIndex > i]
				print(paste("Number of Markers Written into BIG File: ", xx, sep=""))
			}
			if(i >= nmarkers){
				myGeno.backed [(i-maxLine + 1):nmarkers, ] <- -1 * apply(geno[, (i-maxLine + 1):nmarkers], 1, Num.fun) + 3
			}else{
				myGeno.backed [(i-maxLine + 1):i, ] <- -1 * apply(geno[, (i-maxLine + 1):i], 1, Num.fun) + 3
			}
		}
        geno.flush <- flush(myGeno.backed)
        if(!geno.flush){
            stop("flush failed")
        }else{
            print("Preparation for Genotype data is done!")
        }
    }
    
    if(!is.null(fileHMP)){
        #Output Y in same order with Genotype
        #input phenotype
        
        #get the first hmp
        fileHMPCon<-file(description=fileHMP[1], open="r")
        #fileHMPCon<-file(description="hapmap3.hmp.txt", open="r")
        tt<-readLines(fileHMPCon, n=1)
        gg<-readLines(fileHMPCon, n=1, skipNul=1)
        close.connection(fileHMPCon)
        tt2<-unlist(strsplit(tt, sep.hmp))
        gg2<-unlist(strsplit(gg, sep.hmp))
        taxa.g<-as.vector(tt2[-c(1:11)] )
        gg.name<-names(which.max(table(nchar(gg2))))
        bit <- as.numeric(gg.name)
        #length(taxa.g)
        if(is.null(filePhe)){
            taxa=taxa.g
        }else{
            #Output Y in new order
            taxa=intersect(taxa.g, taxa)
            if(length(taxa) == 0){
                print(paste("Phenotype individuals: ", paste(taxa.y[1:5], collapse=", "), "...", sep=""))
                print(paste("Genotype individuals: ", paste(taxa.g[1:5], collapse=", "), "...", sep=""))
                stop("No common individuals between phenotype and genotype!")
            }
            index=match(taxa, myY[, 1], nomatch = 0)
            write.table(myY[index, ], file=paste(out, ".phe", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
            print("Preparation for PHENOTYPE data is done!")
        }

        nFile=length(fileHMP)
        
        #Iteration among file
        print("Output numeric genotype...")
        if(priority == "memory"){
            for (theFile in 1:nFile){
                
                #Open HMP files
                fileHMPCon<-file(description=fileHMP[theFile], open="r")
                
                #Get HMP hearder
                tt<-readLines(fileHMPCon, n=1)
                tt2<-unlist(strsplit(x=tt, split=sep.hmp, fixed=TRUE))
                taxa.g=as.vector(tt2[-c(1:11)] )
                index=match(taxa, taxa.g, nomatch = 0)
                ns=length(index)  #Number of individuals
                
                #Handler maximum number of individuals to output
                if(ns>maxInd)ns=maxInd
                
                #handler of first file
                if(theFile == 1){
                    #Open GD and GM file
                    fileNumCon<-file(description=paste(out, ".numeric", sep=""), open="w")
                    fileMapCon<-file(description=paste(out, ".map", sep=""), open="w")
                    #GM header
                    writeLines("SNP", fileMapCon, sep=sep.num)
                    writeLines("Chrom", fileMapCon, sep=sep.num)
                    writeLines("BP", fileMapCon, sep="\n")
                    #GD header
                    #writeLines(as.character(taxa[1:(ns-1)]), fileNumCon, sep=sep.num)
                    #writeLines(as.character(taxa[ns]), fileNumCon, sep="\n")
                }
                
                #Initialization for iteration within file
                inFile=TRUE
                i=0
                #Iteration within file
                while(inFile){
                    i=i + 1
                    if(i %% 1000 == 0)print(paste("Number of Markers Written into File: ", theFile, ": ", i, sep=""))
                    tt<-readLines(fileHMPCon, n=1)
                    tt2<-unlist(strsplit(x=tt, split=sep.hmp, fixed=TRUE))
                    
                    #Identify end of file
                    if(is.null(tt2[1]))inFile=FALSE
                    if(i>maxRecord)inFile=FALSE
                    
                    if(inFile){
                        #GM
                        rs=tt2[1]
                        chrom=tt2[3]
                        pos=tt2[4]
                        writeLines(as.character(rs), fileMapCon, sep=sep.num)
                        writeLines(as.character(chrom), fileMapCon, sep=sep.num)
                        writeLines(as.character(pos), fileMapCon, sep="\n")
                        
                        #GD
                        GD= MVP.Numericalization(x=tt2[-c(1:11)], bit=bit, effect=SNP.effect, impute=SNP.impute)
                        writeLines(as.character(GD[index[1:(ns-1)]]), fileNumCon, sep=sep.num)
                        writeLines(as.character(GD[index[ns]]), fileNumCon, sep="\n")
                    }#enf of inFIle
                } #end iteration within file
                
                #Close HMP file
                close.connection(fileHMPCon)
            } #end iteration among files
            
            #Close GD and GM file
            close.connection(fileNumCon)
            close.connection(fileMapCon)
        }
        if(priority == "speed"){
            fileNumCon<-file(description=paste(out, ".numeric", sep=""), open="w")
            fileMapCon<-file(description=paste(out, ".map", sep=""), open="w")
            #GM header
            writeLines("SNP", fileMapCon, sep=sep.num)
            writeLines("Chrom", fileMapCon, sep=sep.num)
            writeLines("BP", fileMapCon, sep="\n")
            close.connection(fileNumCon)
            close.connection(fileMapCon)
            for(theFile in 1: nFile){
                fileHMPCon<-file(description=fileHMP[theFile], open="r")
                #Get HMP hearder
                tt<-readLines(fileHMPCon, n=1)
                close.connection(fileHMPCon)
                tt2<-unlist(strsplit(x=tt, split=sep.hmp, fixed=TRUE))
                taxa.g <- as.vector(tt2[-c(1:11)])
                index <- match(taxa, taxa.g, nomatch = 0)
                myFile <- read.delim(fileHMP[theFile], colClasses="character", sep=sep.hmp, head=FALSE, skip=1)
                nM <- nrow(myFile)
                write.table(myFile[, c(1, 3, 4)], paste(out, ".map", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep.num)
                myFile <- myFile[, -c(1:11)];gc()
                myGDx <- apply(myFile[, index], 1, function(x) MVP.Numericalization(x, bit=bit, effect=SNP.effect, impute=SNP.impute))
                myGDx <- t(myGDx)
                write.table(myGDx, paste(out, ".numeric", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep.num)
                rm(myFile);rm(myGDx);gc()
                
                print(paste("File:", fileHMP[theFile], "; Total markers: ", nM, " finished!"))
            }
        }
        print("Preparation for numeric data is done!")
    }else if(!is.null(filePhe)){
        write.table(myY, file=paste(out, ".phe", sep=""), quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = TRUE)
        print("Preparation for PHENOTYPE data is done!")
    }
    
    #map
    if((!is.null(fileMap))|(!is.null(fileHMP))){
    if(is.null(fileMap))	fileMap <- paste(out, ".map", sep="")
        myGM <- read.big.matrix(fileMap, type=type.map, sep=sep.map, head = TRUE)
        #backingfile=paste(out, ".map.bin", sep=""), descriptorfile=paste(out, ".map.desc", sep=""))
        nmarkers <- nrow(myGM)#build file backed matrix for genotype
        #print("Preparation for MAP is done!")
        rm("myGM")
        gc()
        }
    
    #Transfer genotype data to .desc, .bin files
    if((!is.null(fileNum))|(!is.null(fileHMP))){
        print("Output mvp genotype...")
        if(is.null(fileNum)) fileNum <- paste(out, ".numeric", sep="")
        # myGD <- read.big.matrix(fileNum, type=type.geno, sep=sep.num, head = FALSE,
        # backingfile=paste(out, ".geno.bin", sep=""), descriptorfile=paste(out, ".geno.desc", sep=""))
        fileGenoCon <- file(description=fileNum, open="r")
        tt2 <-readLines(fileGenoCon, n=1)
        tt2 <- unlist(strsplit(x=tt2, split=sep.num, fixed=TRUE))
        ns <- length(tt2)
        close.connection(fileGenoCon)
        bck <- paste(out, ".geno.bin", sep="")
        dsc <- paste(out, ".geno.desc", sep="")
        myGeno.backed<-big.matrix(nmarkers, ns, type=type.geno,
        backingfile=bck, descriptorfile=dsc)
        if(priority == "memory"){
            #Initialization for iteration within file
            inGENOFile=TRUE
            i=0
            printN <- unlist(strsplit(x=as.character(nmarkers), split="", fixed=TRUE))
            printIndex <- seq(0, (as.numeric(printN[1]) + 1) * 10^(length(printN)), 1000)[-1]
            #Iteration within file		
            fileGenoCon <- file(description=fileNum, open="r")
            while(inGENOFile){
                i=i + maxLine
                tt<-readLines(fileGenoCon, n=maxLine)
                if(i >= nmarkers){
                    i <- nmarkers
                }
                if(sum(i >= printIndex )!=0){
                    printIndex <- printIndex[printIndex > i]
                    print(paste("Number of Markers Written into MVP File: ", i, sep=""))
                }
                tt<-do.call(rbind, strsplit(x=tt, split=sep.num, fixed=TRUE))
                nn <- nrow(tt)
                #Identify end of file
                if(is.null(tt[1])) inGENOFile=FALSE
                if(inGENOFile){
                    if(i == nmarkers){
                        myGeno.backed [(i-nn + 1):i, ] = tt; rm(tt); gc()
                    }else{
                        myGeno.backed [(i-maxLine + 1):i, ] = tt; rm(tt); gc()
                    }
                }else{
                    geno.flush <- flush(myGeno.backed)
                }
            }
            if(!geno.flush){
                stop("flush failed")
            }else{
                print("Preparation for Genotype File is done!")
            }
            #Close GENO file
            close.connection(fileGenoCon)
        }
        if(priority == "speed"){
            myGeno <- read.big.matrix(fileNum, type=type.geno, head=FALSE, sep=sep.num)
            options(bigmemory.typecast.warning=FALSE)
            myGeno.backed[, ] <- myGeno[, ]
            flush(myGeno.backed)
            rm("myGeno")
        }
        rm(list=c("myGeno.backed"))
        gc()
    }
    
    #Transfer KINSHIP data to .desc, .bin files
    if(!is.logical(fileKin)){
        # myKinship <- read.big.matrix(fileKin, type=type.kin, sep=sep.kin, head = FALSE,
        # backingfile=paste(out, ".kin.bin", sep=""), descriptorfile=paste(out, ".kin.desc", sep=""))
        myKinship <- read.big.matrix(fileKin, head=FALSE, type=type.kin, sep=sep.kin)
        myKinship.backed<-filebacked.big.matrix(nrow(myKinship), ncol(myKinship), type=type.kin,
        backingfile=paste(out, ".kin.bin", sep=""), descriptorfile=paste(out, ".kin.desc", sep=""))
        myKinship.backed[, ] <- myKinship[, ]
        flush(myKinship.backed)
        print("Preparation for KINSHIP is done!")
        rm(list=c("myKinship", "myKinship.backed"))
        #rm(list=c("myKinship"))
        gc()
    }
    #Calculate Kinship matrix
    if((fileKin == TRUE) & (!is.null(fileHMP)|!is.null(fileNum)|!is.null(fileBed))){
        geno.big <- attach.big.matrix(paste(out, ".geno.desc", sep=""))
        print("Calculate KINSHIP using Vanraden method...")
        myKinship <- MVP.K.VanRaden(geno.big, priority=priority, maxLine=maxLine)
        Kin.backed<-big.matrix(nrow(myKinship), ncol(myKinship), type=type.kin, backingfile=paste(out, ".kin.bin", sep=""),
        descriptorfile=paste(out, ".kin.desc", sep=""))
        Kin.backed[, ] <- myKinship[, ]
        flush(Kin.backed)
        print("Preparation for Kinship matrix is done!")
        rm(list=c("Kin.backed", "myKinship"))
        gc()
    }
    
    #Transfer PC data to .desc, .bin files
    if(!is.logical(filePC)){
        myPC <- read.big.matrix(filePC, head=FALSE, type=type.pc, sep=sep.pc)
        PC.backed <- filebacked.big.matrix(nrow(myPC), ncol(myPC), type=type.pc,
        backingfile=paste(out, ".pc.bin", sep=""), descriptorfile=paste(out, ".pc.desc", sep=""))
        PC.backed[ , ] <- myPC[ , ]
        flush(PC.backed)
        print("Preparation for PC matrix is done!")
        rm(list=c("myPC", "PC.backed"))
        gc()
    }
    #Calculate PC matrix
    if((filePC == TRUE) & (!is.null(fileHMP)|!is.null(fileNum)|!is.null(fileBed))){
        geno.big <- attach.big.matrix(paste(out, ".geno.desc", sep=""))
        myPC <- MVP.PCA(geno.big, perc=perc, pcs.keep=pcs.keep)$PCs
        PC.backed<-big.matrix(nrow(myPC), ncol(myPC), type=type.pc, backingfile=paste(out, ".pc.bin", sep=""),
        descriptorfile=paste(out, ".pc.desc", sep=""))
        PC.backed[, ] <- myPC[, ]
        flush(PC.backed)
        print("Preparation for PC matrix is done!")
        rm(list=c("myPC", "PC.backed"))
        gc()
    }

    print("MVP data prepration accomplished successfully!")
}#end of MVP.Data function

`MVP.Numericalization` <-
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
        x[x=="XX"]="N"
        x[x=="--"]="N"
        x[x=="++"]="N"
        x[x=="//"]="N"
        x[x=="NN"]="N"
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
}#end of MVP.Numericalization function
