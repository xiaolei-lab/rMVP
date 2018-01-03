MVP.GLM <-
function(phe, geno, CV=NULL, cpu=2, priority="speed", memo="MVP.GLM", bar=TRUE){
##########################################################################################################
# Object: To perform GWAS with GLM and MLM model and get the P value of SNPs
#
# Input:
# phe: phenotype, n * 2 matrix
# geno: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
# CV: Covariance, design matrix(n * x) for the fixed effects
# cpu: number of cpus used for parallel computation
# memo: a marker on output file name
#
# Output:
# results: a m * 2 matrix, the first column is the SNP effect, the second column is the P values
#
# Authors: Lilin Yin and Xiaolei Liu
# Build date: Aug 30, 2016
# Last update: May 25, 2017
##########################################################################################################
    
    R.ver <- Sys.info()[['sysname']]
    wind <- R.ver == 'Windows'
    taxa <- colnames(phe)[2]
    r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
    math.cpu <- try(getMKLthreads(), silent=TRUE)
    
    n <- ncol(geno)
    m <- nrow(geno)
	
	if(priority=="speed")	geno <- as.matrix(geno)
	
    ys <- as.numeric(as.matrix(phe[,2]))
    
    if(is.null(CV)){
        X0 <- matrix(1, n)
    }else{
        X0 <- cbind(matrix(1, n), CV)
    }
    
    q0 <- ncol(X0)
    iXX <- matrix(NA,q0+1,q0+1)
    Xt <- matrix(NA,n, q0+1)

    y <- matrix(ys)
    X0X0 <- crossprod(X0)

    X0Y <- crossprod(X0,y)

    YY <- crossprod(y)

    X0X0i <- solve(X0X0)

    #parallel function for GLM model
    eff.glm <- function(i){
		if(bar)	print.f(i)
        # if(i%%1000==0){
            # print(paste("****************", i, "****************",sep=""))
        # }
        #if(cpu>1 & r.open){
        #setMKLthreads(math.cpu)
        #}

        SNP <- geno[i, ]
        #Process the edge (marker effects)
            sy <- crossprod(SNP,y)
            ss <- crossprod(SNP)
            xs <- crossprod(X0,SNP)
            
            B21 <- crossprod(xs, X0X0i)
            t2 <- B21 %*% xs
            B22 <- ss - t2
            invB22 <- 1/B22
            NeginvB22B21 <- crossprod(-invB22,B21)
            B21B21 <- crossprod(B21)
            iXX11 <- X0X0i + as.numeric(invB22) * B21B21
            
            #Derive inverse of LHS with partationed matrix
            iXX[1:q0,1:q0] <- iXX11
            iXX[(q0+1),(q0+1)] <- invB22
            iXX[(q0+1),1:q0] <- NeginvB22B21
            iXX[1:q0,(q0+1)] <- NeginvB22B21
            df <- n-q0-1
            rhs <- c(X0Y,sy)
            effect <- crossprod(iXX,rhs)
            ve <- (YY-crossprod(effect,rhs))/df
            effect <- effect[q0+1]
            t.value <- effect/sqrt(iXX[q0+1, q0+1] * ve)
            p <- 2 * pt(abs(t.value), df, lower.tail=FALSE)
        return(list(effect=effect, p=p))
    }
    
    if(cpu == 1){
	    	math.cpu <- try(getMKLthreads(), silent=TRUE)
	    	mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
                try(setMKLthreads(mkl.cpu), silent=TRUE)
		print.f <- function(i){MVP.Bar(i=i, n=m, type="type1", fixed.points=TRUE)}
        results <- lapply(1:m, eff.glm)
	try(setMKLthreads(math.cpu), silent=TRUE)
    }else{
            if(wind){
				print.f <- function(i){MVP.Bar(i=i, n=m, type="type1", fixed.points=TRUE)}
                cl <- makeCluster(getOption("cl.cores", cpu))
                clusterExport(cl, varlist=c("geno", "ys", "X0"), envir=environment())
                Exp.packages <- clusterEvalQ(cl, c(library(bigmemory)))
                results <- parallel::parLapply(cl, 1:m, eff.glm)
                stopCluster(cl)
            }else{
		tmpf.name <- tempfile()
		tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)		
		writeBin(0, tmpf)
		print.f <- function(i){MVP.Bar(n=m, type="type3", tmp.file=tmpf, fixed.points=TRUE)}
                R.ver <- Sys.info()[['sysname']]
                if(R.ver == 'Linux') {
                    math.cpu <- try(getMKLthreads(), silent=TRUE)
                    try(setMKLthreads(1), silent=TRUE)
                }
                results <- parallel::mclapply(1:m, eff.glm, mc.cores=cpu)
		close(tmpf); unlink(tmpf.name); cat('\n');
                if(R.ver == 'Linux') {
                    try(setMKLthreads(math.cpu), silent=TRUE)
                }
            }
    }
    if(is.list(results)) results <- matrix(unlist(results), m, byrow=TRUE)
    #print("****************GLM ACCOMPLISHED****************")
    return(results)
}#end of MVP.GLM function
