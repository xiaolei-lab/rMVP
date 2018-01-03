`MVP.FarmCPU` <- function(phe, geno, map, CV=NULL, priority="speed", P=NULL, method.sub="reward", method.sub.final="reward", method.bin="EMMA", bin.size=c(5e5,5e6,5e7), bin.selection=seq(10,100,10), memo="MVP.FarmCPU", Prior=NULL, ncpus=2, maxLoop=10, threshold.output=.01, converge=1, iteration.output=FALSE, p.threshold=NA, QTN.threshold=NULL, bound=NULL){
    ##############################################################################################
    # Object: Perform GWAS using FarmCPU method
    # Input
    # phe: phenotype
    # geno: genotype, n by m matrix. This is Genotype Data Pure (GD). THERE IS NOT COLUMN FOR TAXA.
    # map: SNP map information
    # CV: covariates
    # priority: speed or memory
    # P: a start p value for each SNP
    # method.sub, method.sub.final: method used in substitution process
    # method.bin: EMMA or FaSTLMM
    # bin.size: window size in genome
    # bin.selection: a vector, how many windows selected
    # memo: a marker on output file name
    # Prior: four columns, SNP name, Chr, Pos, P
    # ncpus: number of cpus used for parallele
    # maxLoop: maximum number of iterations
    # threshold.output: output GWAS results only for SNPs with p value lower than the threshold.output
    # converge: an indicator for when loops stop
    # iteration.output: whether to output results for FarmCPU iterations
    # p.threshold: if all p values in the 1st iteration are bigger than p.threshold, FarmCPU stops
    # QTN.threshold: Only SNPs have a more significant p value than QTN.threshold have chance to be selected as pseudo QTNs
    # bound: maximum number of SNPs selected as pseudo QTNs for each iteration
    # Requirement: Y, GD, and CV should have same taxa order. GD and GM should have the same order on SNPs
    # Output: p values for each marker
    # Authors: Xiaolei Liu and Zhiwu Zhang
    # Date build: Febuary 24, 2013
    # Last update: May 25, 2017
    ##############################################################################################
    #print("--------------------- Welcome to FarmCPU ----------------------------")
	FarmCPU.BIN <-
	function(Y=NULL,GDP=NULL,GM=NULL,CV=NULL,P=NULL,method="EMMA", b=c(5e5,5e6,5e7), s=seq(10,100,10), theLoop=NULL, bound=NULL, ncpus=2){
		#Input: Y - n by 2 matrix with fist column as taxa name and second as trait
		#Input: GDP - n by m+1 matrix. The first colum is taxa name. The rest are m genotype
		#Input: GM - m by 3  matrix for SNP name, chromosome and BP
		#Input: CV - n by t matrix for t covariate variables.
		#Input: P - m by 1 matrix containing probability
		#Input: method - options are "static", "optimum"
		#Input: b - vecter of length>=1 for bin size
		#Input: s - vecter of length>=1 for size of complexity (number of QTNs)
		#Requirement: Y, GDP and CV have same taxa order. GDP and GM have the same order on SNP
		#Requirement: P and GM are in the same order
		#Requirement: No missing data
		#Output: bin - n by s matrix of genotype
		#Output: binmap - s by 3 matrix for map of bin
		#Output: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
		#Relationship: bin=GDP[,c(seqQTN)], binmap=GM[seqQTN,]
		#Authors: Xiaolei Liu and Zhiwu Zhang
		# Last update: March 28, 2017
		##############################################################################
		#print("FarmCPU.BIN Started")
		
		if(is.null(P)) return(list(bin=NULL,binmap=NULL,seqQTN=NULL))
		
		#Set upper bound for bin selection to squareroot of sample size
		n=nrow(Y)
		#bound=round(sqrt(n)/log10(n))
		if(is.null(bound)){
			bound=round(sqrt(n)/sqrt(log10(n)))
		}
		
		s[s>bound]=bound
		s=unique(s[s<=bound]) #keep the within bound
		
		optimumable=(length(b)*length(s)>1)
		if(!optimumable & method=="optimum"){
			method="static"
		}
		
		if(optimumable){
			s[s>bound]=bound
			#print("optimizing possible QTNs...")
			GP=cbind(GM,P,NA,NA,NA)
			mySpecify=MVP.Specify(GI=GM, GP=GP, bin.size=b, inclosure.size=s)
			seqQTN=which(mySpecify$index==TRUE)
		}
		
		#Method of static
		if(method=="static"&optimumable){
			#print("Via static")
			if(theLoop==2){
				b=b[3]
			}else if(theLoop==3){
				b=b[2]
			}else{
				b=b[1]
			}
			s=bound
			s[s>bound]=bound
			print("Optimizing Pseudo QTNs...")
			GP=cbind(GM,P,NA,NA,NA)
			mySpecify=MVP.Specify(GI=GM,GP=GP,bin.size=b,inclosure.size=s)
			seqQTN.save=which(mySpecify$index==TRUE)
		}
		
		#Method of optimum: FaST-LMM
		#============================Optimize by FaST-LMM============================================
		if(method=="FaST-LMM"&optimumable){
			#print("c(bin.size, bin.selection, -2LL, VG, VE)")
			print("Optimizing Pseudo QTNs...")
			count=0
			for (bin in b){
				for (inc in s){
					count=count+1
					GP=cbind(GM,P,NA,NA,NA)
					mySpecify=MVP.Specify(GI=GM,GP=GP,bin.size=bin,inclosure.size=inc)
					seqQTN=which(mySpecify$index==TRUE)
					GK=t(GDP[seqQTN,])
					myBurger=FarmCPU.Burger(Y=Y[,1:2], CV=CV, GK=GK, ncpus=ncpus, method=method)
					myREML=myBurger$REMLs
					myVG=myBurger$vg #it is unused
					myVE=myBurger$ve #it is unused
					print(c(bin,inc,myREML,myVG,myVE))
					#Recoding the optimum GK
					if(count==1){
						seqQTN.save=seqQTN
						LL.save=myREML
						bin.save=bin
						inc.save=inc
						vg.save=myVG  # for genetic variance
						ve.save=myVE  # for residual variance
					}else{
						if(myREML<LL.save){
							seqQTN.save=seqQTN
							LL.save=myREML
							bin.save=bin
							inc.save=inc
							vg.save=myVG  # for genetic variance
							ve.save=myVE  # for residual variance
						}
					} #end of if(count==1)
				}#loop on bin number
			}#loop on bin size
			#seqQTN=seqQTN.save
		}
		
		#Method of optimum: EMMA
		#============================Optimize by EMMA============================================
		if(method=="EMMA"&optimumable){
			#print("c(bin.size, bin.selection, -2LL, VG, VE)")
			print("Optimizing Pseudo QTNs...")
			m <- length(b)*length(s)
			inc.index = rep(c(1:length(s)), length(b))
			
			seqQTN.optimize.parallel <- function(ii){
				bin.index = floor((ii-0.1)/length(s)) + 1
				bin = b[bin.index]
				inc = s[inc.index[ii]]
				GP=cbind(GM,P,NA,NA,NA)
				mySpecify=MVP.Specify(GI=GM,GP=GP,bin.size=bin,inclosure.size=inc)
				seqQTN=which(mySpecify$index==TRUE)
				GK=t(GDP[seqQTN,])
				myBurger=FarmCPU.Burger(Y=Y[,1:2], CV=CV, GK=GK, ncpus=ncpus, method=method)
				myREML=myBurger$REMLs
				myVG=myBurger$vg #it is unused
				myVE=myBurger$ve #it is unused
				print(c(bin,inc,myREML,myVG,myVE))
				return(list(seqQTN=seqQTN,myREML=myREML))
			}
			
			R.ver <- Sys.info()[['sysname']]
			if(R.ver == 'Linux') {
				math.cpu <- try(getMKLthreads(), silent=TRUE)
				try(setMKLthreads(1), silent=TRUE)
			}
			
			llresults <- parallel::mclapply(1:m, seqQTN.optimize.parallel, mc.cores=ncpus)
			
			if(R.ver == 'Linux') {
				try(setMKLthreads(math.cpu), silent=TRUE)
			}
			
			for(i in 1:m){
				if(i == 1){
					seqQTN.save = llresults[[i]]$seqQTN
					myREML.save = llresults[[i]]$myREML
				}else{
					if(llresults[[i]]$myREML < myREML.save){
						seqQTN.save = llresults[[i]]$seqQTN
						myREML.save = llresults[[i]]$myREML
					}
				}
			}
		}
		
		#Method of optimum: GEMMA
		#can not be used to provide REML
		#============================Optimize by EMMA============================================
		if(method=="GEMMA"&optimumable){
			#print("c(bin.size, bin.selection, -2LL, VG, VE)")
			print("Optimizing Pseudo QTNs...")
			m <- length(b)*length(s)
			
			seqQTN.optimize.parallel <- function(ii){
				bin = floor((ii-0.1)/length(s)) + 1
				inc = rep(c(1:length(s)), length(b))
				GP=cbind(GM,P,NA,NA,NA)
				mySpecify=MVP.Specify(GI=GM,GP=GP,bin.size=bin[ii],inclosure.size=inc[ii])
				seqQTN=which(mySpecify$index==TRUE)
				GK=t(GDP[seqQTN,])
				myBurger=FarmCPU.Burger(Y=Y[,1:2], CV=CV, GK=GK, ncpus=ncpus, method=method)
				myREML=myBurger$REMLs
				myVG=myBurger$vg #it is unused
				myVE=myBurger$ve #it is unused
				print(c(bin,inc,myREML,myVG,myVE))
				return(list(seqQTN=seqQTN,myREML=myREML))
			}
			
			llresults <- parallel::mclapply(1:m, seqQTN.optimize.parallel, mc.cores=ncpus)
			
			for(i in 1:m){
				if(i == 1){
					seqQTN.save = llresults[[i]]$seqQTN
					myREML.save = llresults[[i]]$myREML
				}else{
					if(llresults[[i]]$myREML < myREML.save){
						seqQTN.save = llresults[[i]]$seqQTN
						myREML.save = llresults[[i]]$myREML
					}
				}
			}
		}

		return(list(seqQTN=seqQTN.save))
	}#The function FarmCPU.BIN ends here

	MVP.Specify <-
	function(GI=NULL,GP=NULL,bin.size=10000000,inclosure.size=NULL,MaxBP=1e10){
		#Object: To get indicator (TURE or FALSE) for GI based on GP
		#Straitegy
		#1.set bins for all snps in GP
		#2.keep the snp with smallest P value in each bin, record SNP ID
		#3.Search GI for SNP with SNP ID from above
		#4.return the position for SNP selected
		#Input:
		#GI: Data frame with three columns (SNP name, chr and base position)
		#GP: Data frame with seven columns (SNP name, chr and base position, P, MAF,N,effect)
		#Output:
		#theIndex: a vector indicating if the SNPs in GI belong to QTN or not)
		#Authors: Zhiwu Zhang
		#Last update: January 26, 2017
		##############################################################################################
		
		#print("Specification in process...")
		if(is.null(GP))return (list(index=NULL,BP=NULL))
		
		#set inclosure bin in GP
		#Create SNP ID: position+CHR*MaxBP
		ID.GP=as.numeric(as.vector(GP[,3]))+as.numeric(as.vector(GP[,2]))*MaxBP
		
		#Creat bin ID
		bin.GP=floor(ID.GP/bin.size )
		
		#Create a table with bin ID, SNP ID and p value (set 2nd and 3rd NA temporately)
		binP=as.matrix(cbind(bin.GP,NA,NA,ID.GP,as.numeric(as.vector(GP[,4])))  )
		n=nrow(binP)
		
		#Sort the table by p value and then bin ID (e.g. sort p within bin ID)
		binP=binP[order(as.numeric(as.vector(binP[,5]))),]  #sort on P alue
		binP=binP[order(as.numeric(as.vector(binP[,1]))),]  #sort on bin
		
		#set indicator (use 2nd 3rd columns)
		binP[2:n,2]=binP[1:(n-1),1]
		binP[1,2]=0 #set the first
		binP[,3]= binP[,1]-binP[,2]
		
		#Se representives of bins
		ID.GP=binP[binP[,3]>0,]
		
		#Choose the most influencial bins as estimated QTNs
		#Handler of single row
		if(is.null(dim(ID.GP))) ID.GP=matrix(ID.GP,1,length(ID.GP))
		ID.GP=ID.GP[order(as.numeric(as.vector(ID.GP[,5]))),]  #sort on P alue
		
		#Handler of single row (again after reshape)
		if(is.null(dim(ID.GP))) ID.GP=matrix(ID.GP,1,length(ID.GP))
		
		index=!is.na(ID.GP[,4])
		ID.GP=ID.GP[index,4] #must have chr and bp information, keep SNP ID only
		
		if(!is.null(inclosure.size)   ) {
			if(!is.na(inclosure.size)){
				avaiable=min(inclosure.size,length(ID.GP))
				if(avaiable==0){
					ID.GP=-1
				}else{
					ID.GP=ID.GP[1:avaiable] #keep the top ones selected
				}
			}
		}
		
		#create index in GI
		theIndex=NULL
		if(!is.null(GI)){
			ID.GI=as.numeric(as.vector(GI[,3]))+as.numeric(as.vector(GI[,2]))*MaxBP
			theIndex=ID.GI %in% ID.GP
		}
		
		myList=list(index=theIndex,CB=ID.GP)
		return (list(index=theIndex,CB=ID.GP))
	} #end of MVP.Specify

	FarmCPU.LM <-
	function(y,w=NULL,GDP,ncpus=2,npc=0){
		#Object: 1. To quickly sovel LM with one variable substitute multiple times
		#intput: y - dependent variable
		#intput: w - independent variable
		#intput: GDP - independent variable of substitution (GDP)
		#Output: estimate, tvalue, stderr and pvalue ( plus the P value of F test on both A and D)
		#Straitegy: 1. Separate constant covariates (w) and dynamic coveriates (x)
		#Straitegy: 2. Build non-x related only once
		#Straitegy: 3. Use apply to iterate x
		#Straitegy: 4. Derive dominance indicate d from additive indicate (x) mathmaticaly
		#Straitegy: 5. When d is not estimable, continue to test x
		#Authors: Xiaolei Liu and Zhiwu Zhang
		#Start  date: March 1, 2013
		#Last update: March 6, 2013
		##############################################################################################
		#print("FarmCPU.LM started")
		if(is.null(y)) return(NULL)
		if(is.null(GDP)) return(NULL)
		#Constant section (non individual marker specific)
		#---------------------------------------------------------
		#Configration
		nd=20 #number of markes for checking A and D dependency
		N=length(y) #Total number of taxa, including missing ones
		
		if(!is.null(w)){
			w=as.matrix(w,nrow = N)
			nf=ncol(w)
			w=cbind(rep(1,N),w)
			q0=ncol(w)
		}else{
			w=rep(1,N)
			nf=0
			q0=1
		}
		
		print("number of covariates in current loop is:")
		print(nf)
		
		n=N
		if(nd>n)nd=n #handler of samples less than nd
		k=1 #number of genetic effect: 1 and 2 for A and AD respectively
		
		q1=(q0+1) # vecter index for the posistion of genetic effect (a)
		q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
		df=n-q0-k #residual df (this should be varied based on validating d)
		
		iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
		
		ww=crossprodcpp(w)

		wy = crossprod(w,y)
		yy=crossprodcpp(y)
		wwi = solve(ww)
		
		#Statistics on the reduced model without marker
		rhs=wy
		beta <- crossprod(wwi,rhs)
		
		if(npc!=0){
			betapc = beta[2:(npc+1)]
			betapred = beta[-c(1:(npc+1))]
		}else{
			betapc = NULL
			betapred = beta[-1]
		}
		
		m = nrow(GDP)
		
		eff.farmcpu.parallel <- function(i){
			print.f(i)
			# if(i%%1000==0){
				# print(paste("****************", i, "****************",sep=""))
			# }
			x=GDP[i,]
			
			#Process the edge (marker effects)
			xy=crossprod(x,y)
			xx=crossprod(x)
			xw=crossprod(w,x)
			
			B21 <- crossprod(xw, wwi)
			t2 <- B21 %*% xw #I have problem of using crossprod and tcrossprod here
			B22 <- xx - t2
			invB22 = 1/B22
			NeginvB22B21 <- crossprod(-invB22,B21)
			iXX11 <- wwi + as.numeric(invB22) * crossprod(B21)
			
			#Derive inverse of LHS with partationed matrix
			iXX[1:q0,1:q0]=iXX11
			iXX[q1,q1]=invB22
			iXX[q1,1:q0]=NeginvB22B21
			iXX[1:q0,q1]=NeginvB22B21
			
			#statistics
			rhs=c(wy,xy) #the size varied automaticly by A/AD model and validated d
			beta <- crossprod(iXX,rhs)
			df=n-q0-1
			ve=(yy-crossprod(beta,rhs))/df #this is a scaler
			
			#using iXX in the same as above to derive se
			se=sqrt(diag(iXX)*ve)
			tvalue=beta/se
			pvalue <- 2 * pt(abs(tvalue), df,lower.tail = FALSE)
			
			#Handler of dependency between  marker are covariate
			if(!is.na(abs(B22[1,1]))){
				if(abs(B22[1,1])<10e-8)pvalue[]=NA
			}
			B = beta[length(beta)]
			P = pvalue[-1]
			return(list(B=B,P=P))
		}
		print.f <- function(i){MVP.Bar(i=i, n=m, type="type1", fixed.points=TRUE)}
		results <- lapply(1:m, eff.farmcpu.parallel)
		if(is.list(results)) results <- matrix(unlist(results), m, byrow=TRUE)
		return(list(P=results[,-1],betapred=betapred,B=results[,1]))
	} #end of FarmCPU.LM function

	FarmCPU.Burger <-
	function(Y=NULL,CV=NULL,GK=NULL,ncpus=2, method="FaST-LMM"){
		#Object: To calculate likelihood, variances and ratio, revised by Xiaolei based on GAPIT.Burger function from GAPIT package
		#Straitegy: NA
		#Output: P value
		#intput:
		#Y: phenotype with columns of taxa,Y1,Y2...
		#CV: covariate variables with columns of taxa,v1,v2...
		#GK: Genotype data in numerical format, taxa goes to row and snp go to columns. the first column is taxa (same as GAPIT.bread)
		#method: two options for estimating REML, "FaST-LMM" and "EMMA"
		#Authors: Xiaolei Liu and Zhiwu Zhang
		#Last update: Dec 21, 2016
		##############################################################################################
		#print("FarmCPU.Burger in progress...")
		
		if(!is.null(CV)){
			CV=as.matrix(CV)#change CV to a matrix when it is a vector xiaolei changed here
			theCV=as.matrix(cbind(matrix(1,nrow(CV),1),CV)) ###########for FarmCPU
		}else{
			theCV=matrix(1,nrow(Y),1)
		}
		
		#handler of single column GK
		n=nrow(GK)
		m=ncol(GK)
		if(m>2){
			theGK=as.matrix(GK)#GK is pure genotype matrix
		}else{
			theGK=matrix(GK,n,1)
		}
		
		if(method=="FaST-LMM"){
			myFaSTREML=MVP.FaSTLMM.LL(pheno=matrix(Y[,-1],nrow(Y),1), snp.pool=theGK, X0=theCV, ncpus=ncpus)
			REMLs=-2*myFaSTREML$LL
			delta=myFaSTREML$delta
			vg=myFaSTREML$vg
			ve=myFaSTREML$ve
		}
		
		if(method=="EMMA"){
			theGK <- t(theGK)
			K <- MVP.K.VanRaden(M=theGK, priority="speed")
			myEMMAREML <- MVP.EMMA.Vg.Ve(y=matrix(Y[,-1],nrow(Y),1), X=theCV, K=K, cpu=ncpus)
			REMLs=-2*myEMMAREML$REML
			delta=myEMMAREML$delta
			vg=myEMMAREML$vg
			ve=myEMMAREML$ve
		}
		
		#print("FarmCPU.Burger succeed!")
		return (list(REMLs=REMLs,vg=vg,ve=ve,delta=delta))
	} #end of FarmCPU.Burger

	FarmCPU.SUB <-
	function(GM=NULL,GLM=NULL,QTN=NULL,method="mean"){
		#Input: FarmCPU.GLM object
		#Input: QTN - s by 3  matrix for SNP name, chromosome and BP
		#Input: method - options are "penalty", "reward","mean","median",and "onsite"
		#Requirement: P has row name of SNP. s<=t. covariates of QTNs are next to SNP
		#Output: GLM with the last column of P updated by the substituded p values
		#Authors: Xiaolei Liu and Zhiwu Zhang
		# Last update: Febuary 26, 2013
		##############################################################################
		if(is.null(GLM$P)) return(NULL)  #P is required
		if(is.null(QTN)) return(NULL)  #QTN is required
		#print("FarmCPU.SUB Started")
		#print(length(QTN))
		if(length(QTN)==3){
			QTN=QTN[1]
		}else{
			QTN=QTN[,1]
		}
		position=match(QTN, GM[,1], nomatch = 0)
		nqtn=length(position)
		if(is.numeric(GLM$P)){
			GLM$P = as.matrix(GLM$P)
		}
		GLM$B = as.matrix(GLM$B)
		index=(ncol(GLM$P)-nqtn):(ncol(GLM$P)-1)
		spot=ncol(GLM$P)
		if(ncol(GLM$P)!=1){
			if(length(index)>1){
				if(method=="penalty") P.QTN=apply(GLM$P[,index],2,max,na.rm=TRUE)
				if(method=="reward") P.QTN=apply(GLM$P[,index],2,min,na.rm=TRUE)
				if(method=="mean") P.QTN=apply(GLM$P[,index],2,mean,na.rm=TRUE)
				if(method=="median") P.QTN=apply(GLM$P[,index],2,median,na.rm=TRUE)
				if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
			}else{
				if(method=="penalty") P.QTN=max(GLM$P[,index],na.rm=TRUE)
				if(method=="reward") P.QTN=min(GLM$P[,index],na.rm=TRUE)
				if(method=="mean") P.QTN=mean(GLM$P[,index],na.rm=TRUE)
				if(method=="median") P.QTN=median(GLM$P[,index],median,na.rm=TRUE)
				if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
			}
			#replace SNP pvalues with QTN pvalue
			GLM$P[position,spot]=P.QTN
			GLM$B[position,]=GLM$betapred
		}
		return(GLM)
	}#The function FarmCPU.SUB ends here

	FarmCPU.Remove <-
	function(GDP=NULL,GM=NULL,seqQTN=NULL,seqQTN.p=NULL,threshold=.99){
		#Objective: Remove bins that are highly correlated
		#Input: GDP - n by m+1 matrix. The first colum is taxa name. The rest are m genotype
		#Input: GM - m by 3  matrix for SNP name, chromosome and BP
		#Input: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
		#Requirement: GDP and GM have the same order on SNP
		#Output: bin - n by s0 matrix of genotype
		#Output: binmap - s0 by 3 matrix for map of bin
		#Output: seqQTN - s0 by 1 vecter for index of QTN on GM (+1 for GDP column wise)
		#Relationship: bin=GDP[,c(seqQTN)], binmap=GM[seqQTN,], s0<=s
		#Authors: Zhiwu Zhang
		# Last update: March 4, 2013
		##############################################################################
		#print("FarmCPU.Remove Started")
		#print(date())
		
		if(is.null(seqQTN))return(list(bin=NULL,binmap=NULL,seqQTN=NULL))
		seqQTN=seqQTN[order(seqQTN.p)]
		
		hugeNum=10e10
		n=length(seqQTN)
		#fielter bins by physical location
		
		binmap=GM[seqQTN,]
		
		cb=as.numeric(binmap[,2])*hugeNum+as.numeric(binmap[,3])#create ID for chromosome and bp
		cb.unique=unique(cb)
		
		#print("debuge")
		#print(cb)
		#print(cb.unique)
		
		index=match(cb.unique,cb,nomatch = 0)
		seqQTN=seqQTN[index]
		
		#print("Number of bins after chr and bp fillter")
		n=length(seqQTN) #update n
		#print(n)
		#print(date())
		
		#Set sample
		ratio=.1
		maxNum=100000
		m=nrow(GDP) #sample size
		s=ncol(GDP) #number of markers
		
		sampled=s
		if(sampled>maxNum)sampled=maxNum
		
		
		#index=sample(s,sampled)
		index=1:sampled
		
		#This section has problem of turning big.matrix to R matrix
		#It is OK as x is small
		if(is.big.matrix(GDP)){
			x=t(as.matrix(deepcopy(GDP,rows=seqQTN,cols=index) ))
		}else{
			x=t(GDP[seqQTN,index] )
		}
		
		r=cor(as.matrix(x))
		index=abs(r)>threshold
		
		b=r*0
		b[index]=1
		c=1-b
		
		#The above are replaced by following
		c[lower.tri(c)]=1
		diag(c)=1
		bd <- apply(c,2,prod)
		
		position=(bd==1)
		seqQTN=seqQTN[position]
		#============================end of optimum============================================
		seqQTN=seqQTN[!is.na(seqQTN)]
		
		#This section has problem of turning big.matrix to R matrix
		if(is.big.matrix(GDP)){
			bin=t(as.matrix(deepcopy(GDP,rows=seqQTN,) ))
		}else{
			bin=t(GDP[seqQTN,] )
		}
		
		binmap=GM[seqQTN,]
		
		return(list(bin=bin,binmap=binmap,seqQTN=seqQTN))
	}#The function FarmCPU.Remove ends here

	FarmCPU.Prior <-
	function(GM,P=NULL,Prior=NULL){
		#Object: Set prior on existing p value
		#Input: GM - m by 3  matrix for SNP name, chromosome and BP
		#Input: Prior - s by 4  matrix for SNP name, chromosome, BP and Pvalue
		#Input: P - m by 1 matrix containing probability
		#Requirement: P and GM are in the same order, Prior is part of GM except P value
		#Output: P - m by 1 matrix containing probability
		#Authors: Zhiwu Zhang
		# Last update: March 10, 2013
		##############################################################################
		#print("FarmCPU.Prior Started")
		
		if(is.null(Prior)& is.null(P))return(P)
		
		#get prior position
		if(!is.null(Prior)) index=match(Prior[,1],GM[,1],nomatch = 0)
		
		#if(is.null(P)) P=runif(nrow(GM)) #set random p value if not provided (This is not helpful)
		#print("debug set prior  a")
		
		#Get product with prior if provided
		if(!is.null(Prior) & !is.null(P)  )P[index]=P[index]*Prior[,4]
		
		return(P)
	}#The function FarmCPU.Prior ends here
	
    echo=TRUE
    nm=nrow(map)
    if(!is.null(CV)){
        CV=as.matrix(CV)
        npc=ncol(CV)
    }else{
        npc=0
    }
    
    if(is.null(QTN.threshold)){QTN.threshold = 0.01}
    if(!is.na(p.threshold)) QTN.threshold = max(p.threshold, QTN.threshold)
    
    name.of.trait=colnames(phe)[2]
    if(!is.null(memo)) name.of.trait=paste(memo,".",name.of.trait,sep="")
    theLoop=0
    theConverge=0
    seqQTN.save=c(0)
    seqQTN.pre=c(-1)
    isDone=FALSE
    name.of.trait2=name.of.trait
    
    while(!isDone) {
        theLoop=theLoop+1
        print(paste("Current loop: ",theLoop," out of maximum of ", maxLoop, sep=""))
            
        spacer="0"
        if(theLoop>9){
            spacer=""
        }
        if(iteration.output){
            name.of.trait2=paste("Iteration_",spacer,theLoop,".",name.of.trait,sep="")
        }
            
        #Step 2a: Set prior
        myPrior=FarmCPU.Prior(GM=map,P=P,Prior=Prior)

        #Step 2b: Set bins
        if(theLoop<=2){
            myBin=FarmCPU.BIN(Y=phe[,c(1,2)],GD=geno,GM=map,CV=CV,P=myPrior,method=method.bin,b=bin.size,s=bin.selection,theLoop=theLoop,bound=bound,ncpus=ncpus)
        }else{
            myBin=FarmCPU.BIN(Y=phe[,c(1,2)],GD=geno,GM=map,CV=theCV,P=myPrior,method=method.bin,b=bin.size,s=bin.selection,theLoop=theLoop,ncpus=ncpus)
        }
        
        #Step 2c: Remove bin dependency
        #Remove QTNs in LD
        seqQTN=myBin$seqQTN

        if(theLoop==2){
            if(!is.na(p.threshold)){
                if(min(myPrior,na.rm=TRUE)>p.threshold){
                    seqQTN=NULL
                    print("Top snps have little effect, set seqQTN to NULL!")
                    }
                }else{
                    if(min(myPrior,na.rm=TRUE)>0.01/nm){
                        seqQTN=NULL
                        print("Top snps have little effect, set seqQTN to NULL!")
                    }
                }
            }
            
            #when FarmCPU can not work, make a new QQ plot and manhatthan plot
            if(theLoop==2&&is.null(seqQTN)){
                #Report
                P=myGLM$P[,ncol(myGLM$P)]
                P[P==0] <- min(P[P!=0],na.rm=TRUE)*0.01
                results = cbind(myGLM$B, P)
                colnames(results) = c("effect","p")
                break
            }#force to exit for GLM model while seqQTN=NULL and h2=0

            if(!is.null(seqQTN.save)&&theLoop>1){
                if(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN) ) seqQTN=union(seqQTN,seqQTN.save) #Force previous QTNs in the model
            }
            if(theLoop!=1){
                seqQTN.p=myPrior[seqQTN]
                if(theLoop==2){
                    index.p=seqQTN.p<QTN.threshold
                    #if(!is.na(p.threshold)){
                    #index.p=seqQTN.p<p.threshold
                    #}
                    seqQTN.p=seqQTN.p[index.p]
                    seqQTN=seqQTN[index.p]
                    seqQTN.p=seqQTN.p[!is.na(seqQTN)]
                    seqQTN=seqQTN[!is.na(seqQTN)]
                }else{
                    index.p=seqQTN.p<QTN.threshold
                    #if(!is.na(p.threshold)){
                    #index.p=seqQTN.p<p.threshold
                    #}
                    index.p[seqQTN%in%seqQTN.save]=TRUE
                    seqQTN.p=seqQTN.p[index.p]
                    seqQTN=seqQTN[index.p]
                    seqQTN.p=seqQTN.p[!is.na(seqQTN)]
                    seqQTN=seqQTN[!is.na(seqQTN)]
                }
            }

            myRemove=FarmCPU.Remove(GD=geno,GM=map,seqQTN=seqQTN,seqQTN.p=seqQTN.p,threshold=.7)
            
            #Recoding QTNs history
            seqQTN=myRemove$seqQTN
            theConverge=length(intersect(seqQTN,seqQTN.save))/length(union(seqQTN,seqQTN.save))
            circle=(length(union(seqQTN,seqQTN.pre))==length(intersect(seqQTN,seqQTN.pre)))
            
            #handler of initial status
            if(is.null(seqQTN.pre)){circle=FALSE
            }else{
                if(seqQTN.pre[1]==0) circle=FALSE
                if(seqQTN.pre[1]==-1) circle=FALSE
            }

            print("seqQTN")
            print(seqQTN)
            print("scanning...")
            if(theLoop==maxLoop){
                print(paste("Total number of possible QTNs in the model is: ", length(seqQTN),sep=""))
            }
            
            isDone=((theLoop>=maxLoop) | (theConverge>=converge) | circle )
            
            seqQTN.pre=seqQTN.save
            seqQTN.save=seqQTN
            
            #Step 3: Screen with bins
            rm(myBin)
            gc()
            
            theCV=CV
            
            if(!is.null(myRemove$bin)){
                if(length(myRemove$seqQTN) == 1){
                    #myRemove$bin = as.matrix(myRemove$bin)
                    myRemove$bin = t(myRemove$bin)
                }
                theCV=cbind(CV,myRemove$bin)
            }
            myGLM=FarmCPU.LM(y=phe[,2],GDP=geno,w=theCV,ncpus=ncpus,npc=npc)
            
            #Step 4: Background unit substitution
            if(!isDone){
                myGLM=FarmCPU.SUB(GM=map,GLM=myGLM,QTN=map[myRemove$seqQTN,],method=method.sub)
            }else{
                myGLM=FarmCPU.SUB(GM=map,GLM=myGLM,QTN=map[myRemove$seqQTN,],method=method.sub.final)
            }
            P=myGLM$P[,ncol(myGLM$P)]
            P[P==0] <- min(P[P!=0],na.rm=TRUE)*0.01
            results = cbind(myGLM$B, P)
            colnames(results) = c("effect","p")
        } #end of while loop
        #print("****************FarmCPU ACCOMPLISHED****************")
        return(results)
}#The MVP.FarmCPU function ends here

