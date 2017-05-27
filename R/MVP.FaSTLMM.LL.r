`MVP.FaSTLMM.LL` <- function(pheno, snp.pool,X0=NULL,ncpus=2){
    #Evaluation of the maximum likelihood using FaST-LMM method
    #Input: ys, xs, vg, delta, Z, X0, snp.pool
    #pheno: phenotype
    #snp.pool: SNP matrix
    #X0: covariates matrix
    #ncpus: number of CPUs used for parallel computation
    #Output: LL
    #Authors: Qishan Wang, Feng Tian and Zhiwu Zhang
    #Modified by Xiaolei Liu
    #Last update: January 11, 2017
    ################################################################################
    #print("MVP.FaSTLMM.LL started")
    y=pheno
    p=0
    deltaExpStart = -5
    deltaExpEnd = 5
    snp.pool=snp.pool[,]
    if(!is.null(snp.pool)&&var(snp.pool)==0){
        deltaExpStart = 100
        deltaExpEnd = deltaExpStart
    }
    if(is.null(X0)) {
        X0 = matrix(1, nrow(snp.pool), 1)
    }
    X=X0
    #########SVD of X
    K.X.svd <- svd(snp.pool)
    d=K.X.svd$d
    d=d[d>1e-08]
    d=d^2
    U1=K.X.svd$u
    U1=U1[,1:length(d)]
    #handler of single snp
    if(is.null(dim(U1))) U1=matrix(U1,ncol=1)
    n=nrow(U1)
    U1TX=crossprod(U1,X)
    U1TY=crossprod(U1,y)
    yU1TY <- y-U1%*%U1TY
    XU1TX<- X-U1%*%U1TX
    IU = -tcrossprod(U1)
    diag(IU) = rep(1,n) + diag(IU)
    IUX=crossprod(IU,X)
    IUY=crossprod(IU,y)
    #Iteration on the range of delta (-5 to 5 in glog scale)
    delta.range <- seq(deltaExpStart,deltaExpEnd,by=0.1)
    m <- length(delta.range)
    #for (m in seq(deltaExpStart,deltaExpEnd,by=0.1)){
    beta.optimize.parallel <- function(ii){
        #p=p+1
        delta <- exp(delta.range[ii])
        #----------------------------calculate beta-------------------------------------
        #######get beta1
        beta1=0
        for(i in 1:length(d)){
            one=matrix(U1TX[i,], nrow=1)
            beta=crossprod(one,(one/(d[i]+delta)))  #This is not real beta, confusing
            beta1= beta1+beta
        }
        
        #######get beta2
        beta2=0
        for(i in 1:nrow(U1)){
            one=matrix(IUX[i,], nrow=1)
            beta=crossprodcpp(one)
            beta2= beta2+beta
        }
        beta2<-beta2/delta
        
        #######get beta3
        beta3=0
        for(i in 1:length(d)){
            one1=matrix(U1TX[i,], nrow=1)
            one2=matrix(U1TY[i,], nrow=1)
            beta=crossprod(one1,(one2/(d[i]+delta)))
            beta3= beta3+beta
        }
        
        ###########get beta4
        beta4=0
        for(i in 1:nrow(U1)){
            one1=matrix(IUX[i,], nrow=1)
            one2=matrix(IUY[i,], nrow=1)
            beta=crossprod(one1,one2)
            beta4= beta4+beta
        }
        beta4<-beta4/delta
        
        #######get final beta
        zw1 <- geninv(beta1+beta2)
        #zw1 <- try(solve(beta1+beta2))
        #if(inherits(zw1, "try-error")){
        #zw1 <- ginv(beta1+beta2)
        #}
        
        zw2=(beta3+beta4)
        beta=crossprod(zw1,zw2)
        
        #----------------------------calculate LL---------------------------------------
        ####part 1
        part11<-n*log(2*3.14)
        part12<-0
        for(i in 1:length(d)){
            part12_pre=log(d[i]+delta)
            part12= part12+part12_pre
        }
        part13<- (nrow(U1)-length(d))*log(delta)
        part1<- -1/2*(part11+part12+part13)
        
        ######  part2
        part21<-nrow(U1)
        ######part221
        
        part221=0
        for(i in 1:length(d)){
            one1=U1TX[i,]
            one2=U1TY[i,]
            part221_pre=(one2-one1%*%beta)^2/(d[i]+delta)
            part221 = part221+part221_pre
        }
        
        part222=0
        for(i in 1:n){
            one1=XU1TX[i,]
            one2=yU1TY[i,]
            part222_pre=((one2-one1%*%beta)^2)/delta
            part222= part222+part222_pre
        }
        part22<-n*log((1/n)*(part221+part222))
        part2<- -1/2*(part21+part22)
        
        ################# likihood
        LL<-part1+part2
        part1<-0
        part2<-0
        
        return(list(beta=beta,delta=delta,LL=LL))
    }
    #} # end of Iteration on the range of delta (-5 to 5 in glog scale)
    R.ver <- Sys.info()[['sysname']]
    if(R.ver == 'Linux') {
        math.cpu <- try(getMKLthreads(), silent=TRUE)
        try(setMKLthreads(1), silent=TRUE)
    }
    
    llresults <- parallel::mclapply(1:m, beta.optimize.parallel, mc.cores=ncpus)
    
    if(R.ver == 'Linux') {
        try(setMKLthreads(math.cpu), silent=TRUE)
    }
    
    for(i in 1:m){
        if(i == 1){
            beta.save = llresults[[i]]$beta
            delta.save = llresults[[i]]$delta
            LL.save = llresults[[i]]$LL
        }else{
            if(llresults[[i]]$LL > LL.save){
                beta.save = llresults[[i]]$beta
                delta.save = llresults[[i]]$delta
                LL.save = llresults[[i]]$LL
            }
        }
    }
    #--------------------update with the optimum------------------------------------
    beta=beta.save
    delta=delta.save
    LL=LL.save
    
    #--------------------calculating Va and Vem-------------------------------------
    #sigma_a1
    sigma_a1=0
    for(i in 1:length(d)){
        one1=matrix(U1TX[i,], ncol=1)
        one2=matrix(U1TY[i,], nrow=1)
        #sigma_a1_pre=(one2-one1%*%beta)^2/(d[i]+delta)
        sigma_a1_pre=(one2-crossprod(one1,beta))^2/(d[i]+delta)
        sigma_a1= sigma_a1+sigma_a1_pre
    }
    
    ### sigma_a2
    sigma_a2=0
    
    for(i in 1:nrow(U1)){
        one1=matrix(IUX[i,], ncol=1)
        one2=matrix(IUY[i,], nrow=1)
        #sigma_a2_pre<-(one2-one1%*%beta)^2
        sigma_a2_pre<-(one2-crossprod(one1,beta))^2
        sigma_a2= sigma_a2+sigma_a2_pre
    }
    
    sigma_a2<-sigma_a2/delta
    sigma_a<- 1/n*(sigma_a1+sigma_a2)
    sigma_e<-delta*sigma_a
    
    return(list(beta=beta, delta=delta, LL=LL, vg=sigma_a,ve=sigma_e))
}
