MVP.Report <- function(
	MVP,
	col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1", "red"),
	bin.size=1e6,
	bin.max=NULL,
	pch=19,
	band=1,
	cir.band=0.5,
	H=1.5,
	ylim=NULL,
	cex.axis=1,
	plot.type="b",
	multracks=FALSE,
	cex=c(0.5,0.8,1),
	r=0.3,
	xlab="Chromosome",
	ylab=expression(-log[10](italic(p))),
	xaxs="i",
	yaxs="r",
	outward=FALSE,
	threshold = NULL, 
	threshold.col="red",
	threshold.lwd=1,
	threshold.lty=2,
	amplify= TRUE,
	chr.labels=NULL,
	signal.cex = 0.8,
	signal.pch = 19,
	signal.col="red",
	signal.line=NULL,
	cir.chr=TRUE,
	cir.chr.h=0.6,
	cir.chr.col="black",
	cir.legend=TRUE,
	cir.legend.cex=0.8,
	cir.legend.col="grey45",
	LOG10=TRUE,
	box=FALSE,
	conf.int=TRUE,
	conf.int.col="grey",
	file.output=TRUE,
	file="jpg",
	dpi=300
)
{	

	#plot a circle with a radius of r
	circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
	{
		curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
		curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
	}
	
	Densitplot <- function(
		map,
		col=c("darkgreen", "yellow", "red"),
		main="SNP Density",
		bin=1e6,
		band=3,
		width=5,
		legend.len=10,
		legend.max=NULL,
		legend.pt.cex=3,
		legend.cex=1,
		legend.y.intersp=1,
		legend.x.intersp=1
	)
	{
		map <- as.matrix(map)
		map <- map[!is.na(map[, 2]), ]
		map <- map[!is.na(map[, 3]), ]
		map <- map[map[, 2] != 0, ]
		map <- map[map[, 3] != 0, ]
		options(warn = -1)
		max.chr <- max(as.numeric(map[, 2]), na.rm=TRUE)
		if(is.infinite(max.chr))	max.chr <- 0
		map.xy.index <- which(!as.numeric(map[, 2]) %in% c(0 : max.chr))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(map[map.xy.index, 2])
			for(i in 1:length(chr.xy)){
				map[map[, 2] == chr.xy[i], 2] <- max.chr + i
			}
		}
		map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
		chr <- as.numeric(map[, 2])
		pos <- as.numeric(map[, 3])
		chr.num <- unique(chr)
		chorm.maxlen <- max(pos)
		plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
		pos.x <- list()
		col.index <- list()
		maxbin.num <- NULL
		for(i in 1 : length(chr.num)){
			pos.x[[i]] <- pos[which(chr == chr.num[i])]
			cut.len <- round((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
			if(cut.len <= 1){
				col.index[[i]] = 1
			}else{
				cut.r <- cut(pos.x[[i]], cut.len, labels=FALSE)
				eachbin.num <- table(cut.r)
				maxbin.num <- c(maxbin.num, max(eachbin.num))
				col.index[[i]] <- rep(eachbin.num, eachbin.num)
			}
		}
		Maxbin.num <- max(maxbin.num)
		maxbin.num <- Maxbin.num
		if(!is.null(legend.max)){
			maxbin.num <- legend.max
		}
		col=colorRampPalette(col)(maxbin.num)
		for(i in 1 : length(chr.num)){
			polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
				c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
				width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey", border="grey")
			if(!is.null(legend.max)){
				if(legend.max < Maxbin.num){
					col.index[[i]][col.index[[i]] > legend.max] <- legend.max
				}
			}
			segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
			col=col[round(col.index[[i]] * length(col) / maxbin.num)], lwd=1)
		}
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.num[chr.num == max.chr + i] <- chr.xy[i]
			}
		}
		chr.num <- rev(chr.num)
		mtext(at=seq(band, length(chr.num) * band, band),text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=0.6, line=0.2)
		axis(3, at=seq(0, chorm.maxlen, length=10), labels=c(NA, paste(round((seq(0, chorm.maxlen, length=10))[-1] / 1e6, 0), "Mb", sep="")),
			font=1, cex.axis=0.8, tck=0.01, lwd=2, padj=1.2)
		# image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
		# round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
		# t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
		legend.y <- round(seq(0, maxbin.num, length=legend.len))
		len <- legend.y[2]
		legend.y <- seq(0, maxbin.num, len)
		if(!is.null(legend.max)){
			if(legend.max < Maxbin.num){
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, paste(">=", maxbin.num, sep=""))
					legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
				}else{
					legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
					legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
				}
			}else{
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, maxbin.num)
				}
				legend.y.col <- c(legend.y[-1])
			}
		}else{
			if(!maxbin.num %in% legend.y){
				legend.y <- c(legend.y, paste(">", max(legend.y), sep=""))
				legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
			}else{
				legend.y.col <- c(legend.y[-1])
			}
		}
		legend.y.col <- as.numeric(legend.y.col)
		legend(x=(chorm.maxlen + chorm.maxlen/100), y=( -width/2.5 - band * (length(chr.num) - length(chr.num) - 1)), title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=c("grey", col[round(legend.y.col * length(col) / maxbin.num)]),
			cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
	}

	if(class(MVP) == "list"){
		MVP.res <- cbind(MVP$glm.results, MVP$mlm.results, MVP$farmcpu.results)
		Cnames <- colnames(MVP.res)
		Cnames <- Cnames[seq(2, ncol(MVP.res), 2)]
		MVP.res <- MVP.res[, seq(2, ncol(MVP.res), 2)]
		MVP.res[is.na(MVP.res)] <- 1
		MVP.res <- cbind(MVP$map, MVP.res)
		Pmap <- MVP.res
		colnames(Pmap)[-c(1:3)] <- Cnames
	}else{
		Pmap <- MVP
	}
	
	if(sum(plot.type %in% "b")==1) plot.type=c("c","m","q","d")
	
	taxa=colnames(Pmap)[-c(1:3)]
	
	#SNP-Density plot
	if("d" %in% plot.type){
		print("SNP_Density Plotting...")
		if(file.output){
			if(file=="jpg")	jpeg(paste("MVP.SNP_Density.",paste(taxa,collapse="."),".jpg",sep=""), width = 9*dpi,height=7*dpi,res=dpi,quality = 100)
			if(file=="pdf")	pdf(paste("MVP.SNP_Density.",paste(taxa,collapse="."),".pdf",sep=""), width = 9,height=7)
			if(file=="tiff")	tiff(paste("MVP.SNP_Density.",paste(taxa,collapse="."),".tiff",sep=""), width = 9*dpi,height=7*dpi,res=dpi)
		}
		if(!file.output){
			dev.new(width=9, height=7)
		}
		Densitplot(map=Pmap[,c(1:3)], col=col, bin=bin.size, legend.max=bin.max, main=paste("The number of SNPs within ", bin.size/1e6, "Mb window size", sep=""))
		if(file.output)	dev.off()
	}
	if(length(plot.type) !=1 | (!"d" %in% plot.type)){
		#order Pmap by the name of SNP
		#Pmap=Pmap[order(Pmap[,1]),]
		Pmap <- as.matrix(Pmap)

		#delete the column of SNPs names
		Pmap <- Pmap[,-1]
		Pmap <- na.omit(Pmap)

		#scale and adjust the parameters
		cir.chr.h <- cir.chr.h/5
		cir.band <- cir.band/5
		if(!is.null(threshold)){
			threshold.col <- rep(threshold.col,length(threshold))
			threshold.lwd <- rep(threshold.lwd,length(threshold))
			threshold.lty <- rep(threshold.lty,length(threshold))
		}
		if(length(cex)!=3) cex <- rep(cex,3)
		if(!is.null(ylim)){
			if(length(ylim)==1) ylim <- c(0,ylim)
		}

		#get the number of traits
		R=ncol(Pmap)-2

		#replace the non-euchromosome
		options(warn = -1)
		numeric.chr <- as.numeric(Pmap[, 1])
		options(warn = 0)
		max.chr <- max(numeric.chr, na.rm=TRUE)
		if(is.infinite(max.chr))	max.chr <- 0
		map.xy.index <- which(!numeric.chr %in% c(0:max.chr))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(Pmap[map.xy.index, 1])
			for(i in 1:length(chr.xy)){
				Pmap[Pmap[, 1] == chr.xy[i], 1] <- max.chr + i
			}
		}

		Pmap <- matrix(as.numeric(Pmap), nrow(Pmap))

		#order the GWAS results by chromosome and position
		Pmap <- Pmap[order(Pmap[, 1], Pmap[,2]), ]

		#get the index of chromosome
		chr <- unique(Pmap[,1])
		chr.ori <- chr
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
			}
		}

		pvalueT <- as.matrix(Pmap[,-c(1:2)])

		if(LOG10){
			pvalueT[pvalueT <= 0] <- 1
			pvalueT[pvalueT > 1] <- 1
		}

		#set the colors for the plot
		#palette(heat.colors(1024)) #(heatmap)
		#T=floor(1024/max(pvalue))
		#plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
		if(is.vector(col)){
			col <- matrix(col,R,length(col),byrow=T)
		}
		if(is.matrix(col)){
			#try to transform the colors into matrix for all traits
			col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=T)
		}

		Num <- as.numeric(table(Pmap[,1]))
		Nchr <- length(Num)
		N <- NULL

		#set the colors for each traits
		for(i in 1:R){
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			N[i] <- ceiling(Nchr/length(colx))
		}

		#scale the space parameter between chromosomes
		if(!missing(band)){
			band <- floor(band*(dim(pvalueT)[1]/100))
		}else{
			band <- floor(dim(pvalueT)[1]/100)
		}

		#insert the space into chromosomes and return the midpoint of each chromosome
		ticks <- NULL
		NewP <- NULL
		for(j in 1:R){
			pvalue <- pvalueT[,j]
			for(i in 0:(Nchr-1)){
				if (i==0){
					pvalue <- append(pvalue,rep(Inf,band),after=0)
					ticks[i+1] <- band+floor((Num[i+1])/2)
				}else{
					pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
					ticks[i+1] <- band*(i+1)+sum(Num[1:i])+floor((Num[i+1])/2)
				}
			}
			NewP[[j]] <- pvalue
		}

		#merge the pvalues of traits by column
		pvalueT <- do.call(cbind,NewP)
		if(LOG10 == TRUE){
			logpvalueT <- -log10(pvalueT)
		}else{
			pvalueT <- abs(pvalueT)
			logpvalueT <- pvalueT
		}

		Num <- Num+band
		add <- list()
		for(i in 1:R){
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
		}

		TotalN1 <- nrow(pvalueT)
		TotalN <- TotalN1

		signal.line.index <- NULL
		if(!is.null(threshold)){
			if(!is.null(signal.line)){
				for(l in 1:R){
					signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min(threshold)/nrow(pvalueT)))
				}
				signal.line.index <- unique(signal.line.index)
			}
		}
	}
	#plot circle Manhattan
	if("c" %in% plot.type){
		#print("Starting Circular-Manhattan plot!",quote=F)
		if(file.output==TRUE){
			if(file=="jpg")	jpeg(paste("MVP.Circular_Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = 8*dpi,height=8*dpi,res=dpi,quality = 100)
			if(file=="pdf")	pdf(paste("MVP.Circular_Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = 10,height=10)
			if(file=="tiff")	tiff(paste("MVP.Circular_Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = 8*dpi,height=8*dpi,res=dpi)
		}
		if(!file.output){
			dev.new(width=8, height=8)
		}
		par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
		RR <- r+H*R+cir.band*R
		plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
		if(!is.null(signal.line)){
			if(!is.null(signal.line.index)){
				X1chr <- (RR)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
				Y1chr <- (RR)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
				X2chr <- (r)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
				Y2chr <- (r)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
				segments(X1chr,Y1chr,X2chr,Y2chr,lty=2,lwd=signal.line,col="grey")
			}
		}
		for(i in 1:R){
		
			#get the colors for each trait
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			
			#debug
			#print(colx)
			
			print(paste("Circular_Manhattan Plotting ",taxa[i],"...",sep=""))
			pvalue <- pvalueT[,i]
			logpvalue <- logpvalueT[,i]
			if(is.null(ylim)){
				if(LOG10 == TRUE){
					Max <- ceiling(-log10(min(pvalue[pvalue!=0])))
				}else{
					Max <- ceiling(max(pvalue[pvalue!=Inf]))
					if(Max<=1)
					Max <- max(pvalue[pvalue!=Inf])
				}
			}else{
				Max <- ylim[2]
			}
			Cpvalue <- (H*logpvalue/Max)[-c(1:round(band/2))]
			if(outward==TRUE){
				if(cir.chr==TRUE){
					# XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
					circle.plot(myr=RR,lwd=1.5,add=TRUE)
					a <- 0
					
					#plot the boundary which represents the chromosomes

					for(k in 1:length(chr)){
						if(k==1){
							
							#change the axis from right angle into circle format
							X1chr=(RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							Y1chr=(RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							if(is.null(cir.chr.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
							}else{
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
							}
						}else{
							a=a+sum(Pmap[,1]==chr[k-1])
							X1chr=(RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							Y1chr=(RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							if(is.null(cir.chr.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
							}else{
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
							}		
						}
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=H*(-log10(threshold[thr]/max(dim(Pmap))))/Max
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						significantline1=H*(-log10(min(threshold)/max(dim(Pmap))))/Max
					}
				}
				
				#plot the legend for each trait
				if(cir.legend==TRUE){
				
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify == TRUE){
							HX1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							HY1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=rep(rep(colx,N[i]),add[[i]])[which(Cpvalue>=significantline1)])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
							}
						}
					}
				}
				if(cir.chr==TRUE){
					ticks1=1.07*(RR+cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.07*(RR+cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=(0.9*r)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=(0.9*r)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}
			}
			if(outward==FALSE){
				if(cir.chr==TRUE){
					# XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE)
					circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE)

					a=0
					for(k in 1:length(chr)){
						if(k==1){
							X1chr=(2*cir.band+RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							Y1chr=(2*cir.band+RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
							Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}
						}else{
							a=a+sum(Pmap[,1]==chr[k-1])
							X1chr=(2*cir.band+RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							Y1chr=(2*cir.band+RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
							if(is.null(cir.chr.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
							}else{
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
							}	
						}
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=H*(-log10(threshold[thr]/max(dim(Pmap))))/Max
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						significantline1=H*(-log10(min(threshold)/max(dim(Pmap))))/Max
					}
				}
				if(cir.legend==TRUE){
					
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-1)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify == TRUE){
							HX1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							HY1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=rep(rep(colx,N[i]),add[[i]])[which(Cpvalue>=significantline1)])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
							}
							
						}
					}
				}
				if(cir.chr==TRUE){
					ticks1=1.1*(2*cir.band+RR)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.1*(2*cir.band+RR)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						  angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						  text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=1.0*(RR+cir.band)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.0*(RR+cir.band)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						
							#adjust the angle of labels of circle plot
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}	
				}
			}
		}
		if(file.output==TRUE) dev.off()
		#print("Circular-Manhattan has been finished!",quote=F)
	}

	if("m" %in% plot.type){
		if(multracks==FALSE){
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			for(i in 1:R){
				colx=col[i,]
				colx=colx[!is.na(colx)]
				print(paste("Rectangular_Manhattan Plotting ",taxa[i],"...",sep=""))
					if(file.output==TRUE){
						if(file=="jpg")	jpeg(paste("MVP.Rectangular.Manhattan.",taxa[i],".jpg",sep=""), width = 14*dpi,height=5*dpi,res=dpi,quality = 100)
						if(file=="pdf")	pdf(paste("MVP.Rectangular.Manhattan.",taxa[i],".pdf",sep=""), width = 15,height=6)
						if(file=="tiff")	tiff(paste("MVP.Rectangular.Manhattan.",taxa[i],".tiff",sep=""), width = 14*dpi,height=5*dpi,res=dpi)
						par(mar = c(5,6,4,3),xaxs=xaxs,yaxs=yaxs)
					}
					if(file.output==FALSE) {
						dev.new(width = 15, height = 6)
					}
					pvalue=pvalueT[,i]
					logpvalue=logpvalueT[,i]
					if(is.null(ylim)){
						if(!is.null(threshold)){
							if(sum(threshold!=0)==length(threshold)){
								if(LOG10 == TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
								}
							}else{
								if(LOG10 == TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									if(Max<=1)
									#{
										Max=max(max(pvalue[pvalue!=Inf]))
									# }else{
										# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									# }
								}
							}
						}else{
							if(LOG10 == TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								if(Max<=1)
								#{
									Max=max(max(pvalue[pvalue!=Inf]))
								# }else{
									# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								# }
							}
						}
						if(Max<=1){
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
								cex.axis=cex.axis,cex.lab=1.6,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
						}else{
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
								cex.axis=cex.axis,cex.lab=1.6,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
						}
					}else{
						plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
							cex.axis=cex.axis,cex.lab=1.6,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
					}
					if(is.null(chr.labels)){
						axis(1, at=c(0,ticks),cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
					}else{
						axis(1, at=c(0,ticks),cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
					}
					if(is.null(ylim)){
						if(Max>1){
							#print(seq(0,(Max+1),round((Max+1)/10)))
							axis(2,at=seq(0,(Max+1),round((Max+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(Max+1),round((Max+1)/10)))
						}else{
							axis(2,at=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))),cex.axis=cex.axis,font=2,labels=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))))
						}
					}else{
						if(ylim[2]>1){
							axis(2,at=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)))
						}else{
							axis(2,at=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))),cex.axis=cex.axis,font=2,labels=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))))
						}
					}
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							for(thr in 1:length(threshold)){
								h=-log10(threshold[thr]/max(dim(Pmap)))
								# print(h)
								# print(threshold.col[thr])
								# print(threshold.lty[thr])
								# print(threshold.lwd[thr])
								abline(h=h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr])
							}
							if(amplify == TRUE){
								sgline1=-log10(min(threshold)/max(dim(Pmap)))
								HY1=logpvalue[which(logpvalue>=sgline1)]
								HX1=which(logpvalue>=sgline1)
								
								#cover the points that exceed the threshold with the color "white"
								points(HX1,HY1,pch=pch,cex=cex[2],col="white")
								if(is.null(signal.col)){
									points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=rep(rep(colx,N[i]),add[[i]])[HX1])
								}else{
									points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=signal.col)
								}
							}
						}
					}
				if(box==TRUE) box()
				#if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))	abline(v=which(pvalueT[,i] < min(threshold)/max(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
				if(file.output==TRUE)  dev.off()
			}
			#print("Rectangular-Manhattan has been finished!",quote=F)
		}else{
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			#print("Plotting in multiple tracks!",quote=F)
			if(file.output==TRUE){
				if(file=="jpg")	jpeg(paste("Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste("Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = 15,height=6*R)
				if(file=="tiff")	tiff(paste("Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi)
				par(mfcol=c(R,1),mar=c(0, 6+(R-1)*2, 0, 2),oma=c(4,0,4,0),xaxs=xaxs,yaxs=yaxs)
			}
			if(file.output==FALSE){
				dev.new(width = 15, height = 6)
			}
			for(i in 1:R){
				print(paste("Multracks_Rectangular Plotting ",taxa[i],"...",sep=""))
				pvalue=pvalueT[,i]
				logpvalue=logpvalueT[,i]
				if(is.null(ylim)){
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							if(LOG10 ==TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),-log10(min(threshold)/max(dim(Pmap))))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])),-log10(min(threshold)/max(dim(Pmap))))
							}
						}else{
							if(LOG10 == TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								if(Max<=1)
								#{
									Max=max(max(pvalue[pvalue!=Inf]))
								# }else{
									# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								# }
							}	
						}
					}else{
						if(LOG10 == TRUE){
							Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
						}else{
						Max=max(ceiling(max(pvalue[pvalue!=Inf])))
							if(Max<=1)
							#{
								Max=max(max(pvalue[pvalue!=Inf]))
							# }else{
								# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
							# }
						}
					}
					xn <- ifelse(R == 1, R, R * 2/3)
					if(Max<=1){
						plot(logpvalue,pch=pch,cex=cex[2]*xn,col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
							cex.axis=cex.axis*xn,cex.lab=2*xn,font=2,axes=FALSE)
					}else{
						plot(logpvalue,pch=pch,cex=cex[2]*xn,col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
							cex.axis=cex.axis*xn,cex.lab=2*xn,font=2,axes=FALSE)
					}
				}else{
					plot(logpvalue,pch=pch,cex=cex[2]*xn,col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
						cex.axis=cex.axis*xn,cex.lab=2*xn,font=2,axes=FALSE)
				}
				
				#add the names of traits on plot 
				if(Max<=1){
					text(ticks[1],Max+10^(-ceiling(-log10(Max))),labels=taxa[i],adj=0,font=3,cex=xn)
				}else{
					text(ticks[1],Max+1,labels=taxa[i],adj=0,font=3,cex=xn)
				}
				if(i == R){
					if(is.null(chr.labels)){
						axis(1, at=c(0,ticks),cex.axis=cex.axis*xn,font=2,labels=c("Chr",chr.ori),padj=(xn-1)/2)
					}else{
						axis(1, at=c(0,ticks),cex.axis=cex.axis*xn,font=2,labels=c("Chr",chr.labels),padj=(xn-1)/2)
					}
				}
				if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=2,cex=xn)
				if(is.null(ylim)){
					if(Max>1){
						axis(2,at=seq(0,(Max+1),ceiling((Max+1)/10)),cex.axis=cex.axis*xn,font=2,labels=seq(0,(Max+1),ceiling((Max+1)/10)))
					}else{
						axis(2,at=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))),cex.axis=cex.axis*xn,font=2,labels=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))))
					}
				}else{
					if(ylim[2]>1){
						axis(2,at=seq(0,(ylim[2]+1),ceiling((ylim[2]+1)/10)),cex.axis=cex.axis*xn,font=2,labels=seq(0,(ylim[2]+1),ceiling((ylim[2]+1)/10)))
					}else{
						axis(2,at=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))),cex.axis=cex.axis*xn,font=2,labels=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))))
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							h=-log10(threshold[thr]/max(dim(Pmap)))
							abline(h=h,col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
						}
						if(amplify == TRUE){
							sgline1=-log10(min(threshold)/max(dim(Pmap)))
							HY1=logpvalue[which(logpvalue>=sgline1)]
							HX1=which(logpvalue>=sgline1)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=pch,cex=cex[2]*xn,col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2]*xn,col=rep(rep(colx,N[i]),add[[i]])[HX1])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2]*xn,col=signal.col)
							}
						}
					}
				}
				#if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))	abline(v=which(pvalueT[,i] < min(threshold)/max(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
			}
			
			#add the labels of X-axis
			#mtext(xlab,side=1,padj=2.6,font=2,cex=1.6)
			if(file.output==TRUE) dev.off()
			#print("Rectangular-Manhattan has been finished!",quote=F)
		}
	}
		
	if("q" %in% plot.type){
		#print("Starting QQ-plot!",quote=F)
		if(multracks){
			# if(file.output==TRUE){
                # if(file=="jpg")	jpeg(paste("MVP.Multracks.QQplot.",paste(taxa,collapse="."),".jpg",sep=""), width = R*3.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
                # if(file=="pdf")	pdf(paste("MVP.Multracks.QQplot.",paste(taxa,collapse="."),".pdf",sep=""), width = R*3.5,height=3.5)
                # if(file=="tiff")	tiff(paste("MVP.Multracks.QQplot.",paste(taxa,collapse="."),".tiff",sep=""), width = R*3.5*dpi,height=3.5*dpi,res=dpi)
                # par(mfcol=c(1,R),mar = c(0,1,4,1.5),oma=c(3,5,0,0))
			# }else{
				# dev.new(width = 2.5*R, height = 5.5)
			# }
			# for(i in 1:R){
				# print(paste("Multracks_QQ Plotting ",taxa[i],"...",sep=""))		
				# P.values=as.numeric(Pmap[,i+2])
				# P.values=P.values[!is.na(P.values)]
				# if(LOG10==TRUE){
					# P.values=P.values[P.values>0]
					# P.values=P.values[P.values<=1]
					# N=length(P.values)
					# P.values=P.values[order(P.values)]
				# }else{
					# N=length(P.values)
					# P.values=P.values[order(P.values,decreasing=TRUE)]
				# }
				# p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				# log.Quantiles <- -log10(p_value_quantiles)
				# if(LOG10==TRUE){
					# log.P.values <- -log10(P.values)
				# }else{
					# log.P.values <- P.values
				# }
				# # if(i == 1){
					# # plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,max(log.P.values)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = taxa[i])
				# # }else{
					# # plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,max(log.P.values)),xlab =expression(Expected~~-log[10](italic(p))), ylab="", main = taxa[i])
				# # }
				# plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,max(log.P.values)),xlab ="", ylab="", main = taxa[i])
				# #calculate the confidence interval of QQ-plot
				# if(conf.int==TRUE){
					# N1=length(log.Quantiles)
					# c95 <- rep(NA,N1)
					# c05 <- rep(NA,N1)
					# for(j in 1:N1){
						# xi=ceiling((10^-log.Quantiles[j])*N)
						# if(xi==0)xi=1
						# c95[j] <- qbeta(0.95,xi,N-xi+1)
						# c05[j] <- qbeta(0.05,xi,N-xi+1)
					# }
					# index=length(c95):1
					
					# #plot the confidence interval of QQ-plot
					# polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
				# }
				
				# if(!is.null(threshold.col))	abline(a = 0, b = 1, col = threshold.col[1],lwd=2)
				# points(log.Quantiles, log.P.values, col = col[1],pch=19,cex=cex[3])
				# if(!is.null(threshold)){
					# if(sum(threshold!=0)==length(threshold)){
						# thre.line=-log10(min(threshold)/N)
						# if(amplify==TRUE){
							# thre.index=which(log.P.values>=thre.line)
							# if(length(thre.index)!=0){
							
								# #cover the points that exceed the threshold with the color "white"
								# points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
								# if(is.null(signal.col)){
									# points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch,cex=signal.cex)
								# }else{
									# points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col,pch=signal.pch,cex=signal.cex)
								# }
							# }
						# }
					# }
				# }
			# }
			# if(file.output==TRUE) dev.off()
				if(file.output==TRUE){
					if(file=="jpg")	jpeg(paste("MVP.Multraits.QQplot.jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
					if(file=="pdf")	pdf(paste("MVP.Multraits.QQplot.pdf",sep=""), width = 5.5,height=5.5)
					if(file=="tiff")	tiff(paste("MVP.Multraits.QQplot.tiff",sep=""), width = 5.5*dpi,height=5.5*dpi, res=dpi)
					par(mar = c(5,5,4,2))
				}else{
					dev.new(width = 5.5, height = 5.5)
				}
				p_value_quantiles=(1:nrow(Pmap))/(nrow(Pmap)+1)
				log.Quantiles <- -log10(p_value_quantiles)
                Pmap.min <- Pmap[,3:(R+2)]
				YlimMax <- -log10(min(Pmap.min[Pmap.min > 0])); gc()
                #print("YlimMax")
                #print(YlimMax)
				plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0, YlimMax),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = "QQplot")
				for(i in 1:R){
					print(paste("Multraits.QQ Plotting ",taxa[i],"...",sep=""))
					P.values=as.numeric(Pmap[,i+2])
					P.values=P.values[!is.na(P.values)]
					if(LOG10==TRUE){
						P.values=P.values[P.values>=0]
						P.values=P.values[P.values<=1]
						N=length(P.values)
						P.values=P.values[order(P.values)]
					}else{
						N=length(P.values)
						P.values=P.values[order(P.values,decreasing=TRUE)]
					}
					if(LOG10==TRUE){
						log.P.values <- -log10(P.values)
					}else{
						log.P.values <- P.values
					}
							
					# calculate the confidence interval of QQ-plot
					if((i == 1) & conf.int==TRUE){
						N1=length(log.Quantiles)
						c95 <- rep(NA,N1)
						c05 <- rep(NA,N1)
						for(j in 1:N1){
							xi=ceiling((10^-log.Quantiles[j])*N)
							if(xi==0)xi=1
							c95[j] <- qbeta(0.95,xi,N-xi+1)
							c05[j] <- qbeta(0.05,xi,N-xi+1)
						}
						index=length(c95):1
							
						# plot the confidence interval of QQ-plot
						polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
					}
						
					if((i == 1) & !is.null(threshold.col))	abline(a = 0, b = 1, col = threshold.col[1],lwd=2)
					points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
						
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							thre.line=-log10(min(threshold)/N)
							if(amplify==TRUE){
								thre.index=which(log.P.values>=thre.line)
								if(length(thre.index)!=0){
								
									# cover the points that exceed the threshold with the color "white"
									points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
									if(is.null(signal.col)){
										points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch,cex=signal.cex)
									}else{
										points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col,pch=signal.pch,cex=signal.cex)
									}
								}
							}
						}
					}
				}
				legend("topleft",taxa,col=t(col)[1:R],lwd=2,text.font=6)
				if(file.output==TRUE) dev.off()
		}else{
			for(i in 1:R){
				print(paste("Q_Q Plotting ",taxa[i],"...",sep=""))
				if(file.output==TRUE){
					if(file=="jpg")	jpeg(paste("MVP.QQplot.",taxa[i],".jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
					if(file=="pdf")	pdf(paste("MVP.QQplot.",taxa[i],".pdf",sep=""), width = 5.5,height=5.5)
					if(file=="tiff")	tiff(paste("MVP.QQplot.",taxa[i],".tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
					par(mar = c(5,5,4,2))
				}else{
					dev.new(width = 5.5, height = 5.5)
				}
				P.values=as.numeric(Pmap[,i+2])
				P.values=P.values[!is.na(P.values)]
				if(LOG10==TRUE){
					P.values=P.values[P.values>0]
					P.values=P.values[P.values<=1]
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				if(LOG10==TRUE){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}
				plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,max(log.P.values)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))
				
				#calculate the confidence interval of QQ-plot
				if(conf.int==TRUE){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N)
						if(xi==0)xi=1
						c95[j] <- qbeta(0.95,xi,N-xi+1)
						c05[j] <- qbeta(0.05,xi,N-xi+1)
					}
					index=length(c95):1
					
					#plot the confidence interval of QQ-plot
					polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
				}
				
				if(!is.null(threshold.col))	abline(a = 0, b = 1, col = threshold.col[1],lwd=2)
				points(log.Quantiles, log.P.values, col = col[1],pch=19,cex=cex[3])
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						thre.line=-log10(min(threshold)/N)
						if(amplify==TRUE){
							thre.index=which(log.P.values>=thre.line)
							if(length(thre.index)!=0){
							
								#cover the points that exceed the threshold with the color "white"
								points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
								if(is.null(signal.col)){
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch,cex=signal.cex)
								}else{
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col,pch=signal.pch,cex=signal.cex)
								}
							}
						}
					}
				}
				if(file.output==TRUE) dev.off()
			}
		}
	}
}

MVP.Hist <- function(
phe,
col=c("dodgerblue4","olivedrab4","violetred","darkgoldenrod1","purple4"),
breakNum=15,
file="pdf",
dpi=300
)
{
		if(file=="jpg")	jpeg(paste("MVP.Phe_Distribution.",paste(colnames(phe)[2],collapse="."),".jpg",sep=""), width = 6*dpi,height=6*dpi,res=dpi,quality = 100)
		if(file=="pdf")	pdf(paste("MVP.Phe_Distribution.",paste(colnames(phe)[2],collapse="."),".pdf",sep=""), width = 6,height=6)
		if(file=="tiff")	tiff(paste("MVP.Phe_Distribution.",paste(colnames(phe)[2],collapse="."),".tiff",sep=""), width = 6*dpi,height=6*dpi,res=dpi)
		Breaks <- seq(min(phe[, 2]), max(phe[, 2]), length=breakNum)
		xx <- hist(phe[, 2], breaks=Breaks,xlab="",ylab="Density", freq=FALSE, col=colorRampPalette(col)(breakNum), font=2, font.lab=2, main=paste("Distribution of ", colnames(phe)[2], sep=""))
		maxY <- max(max(xx$density),  max(density(phe[, 2])$y))
		hist(phe[, 2], breaks=Breaks,xlab="",ylab="Density", ylim=c(0, maxY), freq=FALSE, col=colorRampPalette(col)(breakNum), font=2, font.lab=2, main=paste("Distribution of ", colnames(phe)[2], sep=""))
	    lines(density(phe[, 2]), lwd=2)
	
		if(length(phe[, 2]) <= 5000){
			norm.p <- round(shapiro.test(phe[, 2])$p, 4)
			test.method <- "Shapiro-Wilk"
		}else{
			norm.p <- round(ks.test(phe[, 2],"pnorm")$p, 4)
			test.method <- "Kolmogorov-Smirnov"
		}
		text(xx$breaks[1], y=maxY*0.95, labels=paste("Mean: ", round(mean(phe[, 2]), 2), sep=""), font=2, adj=0)
		text(xx$breaks[1], y=maxY*0.9, labels=paste("Sd: ", round(sd(phe[, 2]), 2), sep=""), font=2, adj=0)
		text(xx$breaks[1], y=maxY*0.85, labels=paste(test.method, ": ", norm.p, sep=""), font=2, adj=0)
	    dev.off()
}

MVP.PCAplot <- function(
PCA,
col=c("darkgreen","darkmagenta","orange"),
pch=c(16,17,8),
Ncluster=3,
plot3D=TRUE,
file="pdf",
dpi=300
)
{
	if(!is.null(col)){
		col <- rep(col, Ncluster)
	}
	if(!is.null(pch)){
		pch <- rep(pch, Ncluster)
	}		
	print("PCA plot2d...")
	if(file=="jpg")	jpeg("MVP.PCA_2D.jpg", width = 6*dpi,height=6*dpi,res=dpi,quality = 100)
	if(file=="pdf")	pdf("MVP.PCA_2D.pdf", width = 6,height=6)
	if(file=="tiff")	tiff("MVP.PCA_2D.tiff", width = 6*dpi,height=6*dpi,res=dpi)			
	kc <- kmeans(PCA[,c(1,2)], Ncluster)
	plot(PCA[,1],PCA[,2],pch=pch[kc$cluster],col=col[kc$cluster],font=2,font.lab=2,xlab="PC1",ylab="PC2")
    #kc <- kmeans(PCA[,c(2,3)], Ncluster)
    #plot(PCA[,2],PCA[,3],pch=pch[kc$cluster],col=col[kc$cluster],font=2,font.lab=2,xlab="PC2",ylab="PC3")
    #kc <- kmeans(PCA[,c(1,3)], Ncluster)
    #plot(PCA[,1],PCA[,3],pch=pch[kc$cluster],col=col[kc$cluster],font=2,font.lab=2,xlab="PC1",ylab="PC3")
	dev.off()
	if(plot3D){
		print("PCA plot3d...")
		par3d(cex=0.8,windowRect=c(100,100,600,600),font=4,userMatrix=matrix(c(0.88,0.47,-0.07,0,-0.14,0.40,0.90,0,0.45,-0.78,0.42,0,0,0,0,1),4,byrow=TRUE))
		kc <- kmeans(PCA, Ncluster)
		plot3d(PCA,col=col[kc$cluster],type="s",size=1.5,top=FALSE)
		#rgl.postscript("MVP.PCA_3D.pdf","pdf")
		dir.create('PCA plot3D')
		setwd('./PCA plot3D')
		M <- par3d("userMatrix")
		movie3d(par3dinterp(time=(0:2)*0.75,userMatrix=list(M,rotate3d(M, pi, 1, 0, 0), rotate3d(M, pi, 0, 1, 0))),fps=30,duration=3,movie="MVP.3D",dir=getwd(),top=FALSE)
		movie3d(spin3d(rpm=20),duration=3,fps=15,frames="MVP.3D",convert=FALSE,dir=getwd(),top=FALSE)
		setwd("../")
	}
}

MVP.Bar <- 
function(
	i, n, type=c("type1", "type2", "type3"), symbol="-", tmp.file=NULL, symbol.head="|", symbol.tail=">" ,fixed.points=TRUE, points=seq(0,100,1), symbol.len=70
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: show the rate of progress for loop 	 														 #
#	 																									 #
# Input:	 																							 #
# i: the current loop number	 																		 #
# n: max loop number	 																				 #
# type: type1 for "for" function, type2 for "mclapply" function											 #
# tmp.file: the opened file of "fifo" function															 #
# symbol: the symbol for the rate of progress	 														 #
# symbol.head: the head for the bar																		 #
# symbol.tail: the tail for the bar																		 #
# fixed.points: whether use the setted points which will be printed	 									 #
# points: the setted points which will be printed	 													 #
# symbol.len: the total length of progress bar		 													 #
#--------------------------------------------------------------------------------------------------------#
	switch(
		match.arg(type), 
		"type1"={
			if(fixed.points){
				point.index <- points
				point.index <- point.index[point.index > floor(100*(i-1)/n)]
				if(floor(100*i/n) %in% point.index){
					if(floor(100*i/n) != max(point.index)){
						print.len <- floor(symbol.len*i/n)
						cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							paste(rep(" ", symbol.len-print.len), collapse=""),
							sprintf("%.2f%%", 100*i/n), sep="")
						)
					}else{
						print.len <- floor(symbol.len*i/n)
						cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							sprintf("%.2f%%", 100*i/n), "\n", sep="")
						)
					}
				}
			}else{
				if(i < n){
					print.len <- floor(symbol.len*i/n)
					cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						paste(rep(" ", symbol.len-print.len), collapse=""),
						sprintf("%.2f%%", 100*i/n), sep="")
					)
				}else{
					print.len <- floor(symbol.len*i/n)
					cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						sprintf("%.2f%%", 100*i/n), "\n", sep="")
					)
				}
			}
		},
		"type2"={
			if(inherits(parallel:::mcfork(), "masterProcess")) {
				progress <- 0.0
				while(progress < n && !isIncomplete(tmp.file)){
					msg <- readBin(tmp.file, "double")
					progress <- progress + as.numeric(msg)
					print.len <- round(symbol.len * progress / n)
					if(fixed.points){
						if(progress %in% round(points * n / 100)){
							cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							paste(rep(" ", symbol.len-print.len), collapse=""),
							sprintf("%.2f%%", progress * 100 / n), sep=""))
						}
					}else{
						cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						paste(rep(" ", symbol.len-print.len), collapse=""),
						sprintf("%.2f%%", progress * 100 / n), sep=""))
					}
				}
				parallel:::mcexit()
			}
		},
		"type3"={
			progress <- readBin(tmp.file, "double") + 1
			writeBin(progress, tmp.file)
			print.len <- round(symbol.len * progress / n)
			if(fixed.points){
				if(progress %in% round(points * n / 100)){
					cat(paste("\r", 
					paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
					paste(rep(" ", symbol.len-print.len), collapse=""),
					sprintf("%.2f%%", progress * 100 / n), sep=""))
				}
			}else{
				cat(paste("\r", 
				paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
				paste(rep(" ", symbol.len-print.len), collapse=""),
				sprintf("%.2f%%", progress * 100 / n), sep=""))
			}
		}
	)
}

