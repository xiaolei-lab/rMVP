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


#' MVP.Report
#' 
#' @param MVP a dataframe or list, at least four columns. The first column is the name of SNP, the second column is the chromosome of SNP, the third column is the position of SNP, and the remaining columns are the P-value of each trait(Note:each trait a column). 
#' @param col a vector or a matrix, if "col" is a vector, each circle use the same colors, it means that the same chromosome is drewed in the same color, the colors are not fixed, one, two, three or more colors can be used, if the length of the "col" is shorter than the length the chromosome, then colors will be applied circularly. If "col" is a matrix, the row is the number of circles(traits), the columns are the colors that users want to use for different circles, each circle can be plotted in different number of colors, the missing value can be replaced by NA. For example: col=matrix(c("grey30","grey60",NA,"red","blue","green","orange",NA,NA),3,3,byrow=T). 
#' @param bin.size the size of bin for SNP_density plot. 
#' @param bin.range a vector, c(min, max). The min/max value of legend of SNP_density plot, the bin whose SNP number is smaller/bigger than 'bin.range' will be use the same color. 
#' @param pch a number, the type for the points or for traits of multi-traits Manhattan plot, is the same with "pch" in <plot>. 
#' @param band a number, the space between chromosomes, the default is 1(if the band equals to 0, then there would be no space between chromosomes). 
#' @param H a number, the height for each circle, each circle represents a trait, the default is 1. 
#' @param ylim a vector, the range of Y-axis when plotting the two type of Manhattan plots, is the same with "ylim" in <plot>. 
#' @param cex.axis a number, controls the size of ticks' numbers of X/Y-axis and the size of labels of circle plot. 
#' @param lwd.axis a number, controls the width of X/Y-axis lines. 
#' @param cex.lab a number, controls the size of labels of X/Y-axis. 
#' @param plot.type a character or vector, only "d", "c", "m", "q" or "b" can be used. if plot.type="d", SNP density will be plotted; if plot.type="c", only circle-Manhattan plot will be plotted; if plot.type="m",only Manhattan plot will be plotted; if plot.type="q",only Q-Q plot will be plotted;if plot.type="b", both circle-Manhattan, Manhattan and Q-Q plots will be plotted; if plot.type=c("m","q"), Both Manhattan and Q-Q plots will be plotted. 
#' @param multracks a logical,if multracks=FALSE, all Manhattan plots will be drew in separated files, if it is TRUE, all Manhattan plots will be plotted in only one file. 
#' @param cex a number or a vector, the size for the points, is the same with "size" in <plot>, and if it is a vector, the first number controls the size of points in circle plot(the default is 0.5), the second number controls the size of points in Manhattan plot(the default is 1), the third number controls the size of points in Q-Q plot(the default is 1) 
#' @param r a number, the radius for the circle(the inside radius), the default is 1. 
#' @param xlab a character, the labels for x axis. 
#' @param ylab a character, the labels for y axis. 
#' @param xaxs a character, The style of axis interval calculation to be used for the x-axis. Possible values are "r", "i", "e", "s", "d". The styles are generally controlled by the range of data or xlim, if given. 
#' @param yaxs a character, The style of axis interval calculation to be used for the y-axis. See xaxs above.. 
#' @param outward logical, if outward=TRUE,then all points will be plotted from inside to outside for circular Manhattan plot. 
#' @param threshold a number or vector, the significant threshold. For example, Bonfferoni adjustment method: threshold=0.01/nrow(Pmap). More than one significant line can be added on the plots, if threshold=0 or NULL, then the threshold line will not be added. 
#' @param threshold.col a character or vector, the colour for the line of threshold levels. 
#' @param threshold.lwd a number or vector, the width for the line of threshold levels. 
#' @param threshold.lty a number or vector, the type for the line of threshold levels. 
#' @param amplify logical, CMplot can amplify the significant points, if amplify=T, then the points bigger than the minimal significant level will be amplified, the default: amplify=TRUE. 
#' @param signal.cex a number, if amplify=TRUE, users can set the size of significant points. 
#' @param signal.pch a number, if amplify=TRUE, users can set the shape of significant points. 
#' @param signal.col a character, if amplify=TRUE, users can set the colour of significant points, if signal.col=NULL, then the colors of significant points will not be changed. 
#' @param signal.line a number, the width of the lines of significant SNPs cross the circle. 
#' @param highlight a vector, names of SNPs which need to be highlighted. 
#' @param highlight.cex a number or vector, the size of points for SNPs which need to be highlighted. 
#' @param highlight.pch a number or vector, the pch of points for SNPs which need to be highlighted. 
#' @param highlight.col a number or vector, the col of points for SNPs which need to be highlighted. 
#' @param chr.labels a vector, the labels for the chromosomes of density plot and circle-Manhattan plot. 
#' @param chr.den.col a character or vector or NULL, the colour for the SNP density. If the length of parameter 'chr.den.col' is bigger than 1, SNP density that counts the number of SNP within given size('bin.size') will be plotted around the circle. If chr.den.col=NULL, the density bar will not be attached on the bottom of manhattan plot.
#' @param cir.band a number, the space between circles, the default is 1.
#' @param cir.chr logical, a boundary that represents chromosomes will be plotted on the periphery of a circle, the default is TRUE.
#' @param cir.chr.h a number, the width for the boundary, if cir.chr=FALSE, then this parameter will be useless.
#' @param cir.legend logical, whether to add the legend of each circle.
#' @param cir.legend.cex a number, the size of the number of legend.
#' @param cir.legend.col a character, the color of the axis of legend.
#' @param LOG10 logical, whether to change the p-value into log10(p-value).
#' @param box logical, this function draws a box around the current Manhattan plot.
#' @param conf.int logical, whether to plot confidence interval on QQ-plot.
#' @param file.output a logical, users can choose whether to output the plot results.
#' @param outpath Only when file.output = TRUE, determines the path of the output file
#' @param file.type a character, users can choose the different output formats of plot, so for, "jpg", "pdf", "tiff" can be selected by users.
#' @param dpi a number, the picture resolution for .jpg and .tiff files. The default is 300.
#' @param height the height of output files.
#' @param width the width of output files.
#' @param memo add a character to the output file name.
#' @param verbose whether to print the reminder. 
#'
#' @export
#' @return Output files
#' @examples
#' data(pig60K, package = "rMVP")
#' 
#' MVP.Report(pig60K[,c(1:3, 5)], plot.type="m",
#'     threshold=0.05/nrow(pig60K), file.output=FALSE)
#' 

MVP.Report <- function(
	MVP,
	col=c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
	bin.size=1e6,
	bin.range=NULL,
	pch=19,
	band=1,
	H=1.5,
	ylim=NULL,
	cex.axis=1,
	lwd.axis=1.5,
	cex.lab=1.5,
	plot.type="b",
	multracks=FALSE,
	cex=c(0.5,1,1),
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
	amplify= FALSE,
	signal.cex = 1.5,
	signal.pch = 19,
	signal.col="red",
	signal.line=1,
	highlight=NULL,
	highlight.cex=1.5,
	highlight.pch=19,
	highlight.col="green",
	chr.labels=NULL,
	chr.den.col="black",
	cir.band=1,
	cir.chr=TRUE,
	cir.chr.h=1.5,
	cir.legend=TRUE,
	cir.legend.cex=0.6,
	cir.legend.col="black",
	LOG10=TRUE,
	box=FALSE,
	conf.int=TRUE,
	file.output=TRUE,
	outpath=getwd(),
	file.type="jpg",
	dpi=300,
	height=NULL,
	width=NULL,
	memo="",
	verbose=TRUE
)
{	

	#plot a circle with a radius of r
	circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
	{
		curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
		curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
	}
	
	max_ylim <- function(x){
		if(x == 0) return(x)
		if(abs(x) >= 1){
			return(ceiling(x))
		}else{
			if(x < 0){
				digit <- 10^(-ceiling(-log10(abs(x))))
				return(-(floor(abs(x) / digit - 1) * digit))
			}else{
				digit <- 10^(-ceiling(-log10(x)))
				return((floor(x / digit + 1) * digit))
			}
		}
	}

	min_ylim <- function(x){
		if(x == 0) return(x)
		if(abs(x) >= 1){
			return(floor(x))
		}else{
			if(x < 0){
				digit <- 10^(-ceiling(-log10(abs(x))))
				return(-(floor(abs(x) / digit + 1) * digit))
			}else{
				digit <- 10^(-ceiling(-log10(x)))
				return((floor(x / digit - 1) * digit))
			}
		}
	}

	min_no_na <- function(x){
		return(min(x, na.rm=TRUE))
	}

	max_no_na <- function(x){
		return(max(x, na.rm=TRUE))
	}

	Densitplot <- function(
		map,
		col=c("darkgreen", "yellow", "red"),
		main="SNP Density",
		chr.labels=NULL, 
		bin=1e6,
		band=3,
		width=5,
		legend.len=10,
		legend.max=NULL,
		legend.min=1,
		legend.pt.cex=3,
		legend.cex=1,
		legend.y.intersp=1,
		legend.x.intersp=1,
		plot=TRUE
	)
	{
		if(is.null(legend.min))	legend.min = 1
		if(is.null(col) | length(col) == 1){col=c("darkgreen", "yellow", "red")}
		map <- as.matrix(map)
		map <- map[!is.na(map[, 2]), ]
		map <- map[!is.na(as.numeric(map[, 3])), ]
		map <- map[map[, 2] != 0, ]
		#map <- map[map[, 3] != 0, ]
		suppressWarnings(max.chr <- max(as.numeric(map[, 2]), na.rm=TRUE))
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
		if(plot)	plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
		pos.x <- list()
		col.index <- list()
		maxbin.num <- NULL
		for(i in 1 : length(chr.num)){
			pos.x[[i]] <- pos[which(chr == chr.num[i])]
			cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
			if(cut.len <= 1){
				maxbin.num <- c(maxbin.num,length(pos.x[[i]]))
               	col.index[[i]] <- rep(length(pos.x[[i]]), length(pos.x[[i]]))
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
		col=colorRampPalette(col)(maxbin.num - legend.min + 1)
		col.seg=NULL
		for(i in 1 : length(chr.num)){
			if(plot)	polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
				c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
				width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey", border="grey")
			if(!is.null(legend.max)){
				if(legend.max < Maxbin.num){
					col.index[[i]][col.index[[i]] > legend.max] <- legend.max
				}
			}
			col.index[[i]][col.index[[i]] < legend.min] <- legend.min
			col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
			if(plot)	segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
			col=col[round(col.index[[i]] * length(col) / (maxbin.num - legend.min + 1)) - legend.min + 1], lwd=1)
		}
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.num[chr.num == max.chr + i] <- chr.xy[i]
			}
		}
		chr.num <- rev(chr.num)
		if(plot){
			if(!is.null(chr.labels)){
				mtext(at=seq(band, length(chr.num) * band, band),text=chr.labels, side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2)
			}else{
				if(max.chr == 0)	mtext(at=seq(band, length(chr.num) * band, band),text=chr.num, side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2)
				if(max.chr != 0)	mtext(at=seq(band, length(chr.num) * band, band),text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2)
			}
		}
		if(plot)	axis(3, at=seq(0, chorm.maxlen, length=10), labels=paste(round((seq(0, chorm.maxlen, length=10)) / 1e6, 0), "Mb", sep=""),
			font=1, cex.axis=cex.axis*0.8, tck=0.01, lwd=2, padj=1.2)
		# image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
		# round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
		# t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
				
		if(maxbin.num <= legend.len)	legend.len <- maxbin.num		
		
		legend.y <- round(seq(legend.min, maxbin.num, length=legend.len))
		legend.y <- unique(legend.y)
		len <- ifelse(length(legend.y)==1, 1, legend.y[2] - legend.y[1])
		legend.y <- seq(legend.min, maxbin.num, len)
		if(!is.null(legend.max)){
			if(legend.max < Maxbin.num){
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, paste(">=", maxbin.num, sep=""))
					legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
				}else{
					legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
					legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
				}
			}else{
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, maxbin.num)
				}
				legend.y.col <- c(legend.y)
			}
		}else{
			if(!maxbin.num %in% legend.y){
				legend.y <- c(legend.y, paste(">", max(legend.y), sep=""))
				legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
			}else{
				legend.y.col <- c(legend.y)
			}
		}
		if(legend.min != 1){
			legend.y[1] <- paste("<=", legend.min, sep="")
		}
		legend.y <- c(0, legend.y)
		legend.y.col <- as.numeric(legend.y.col)
		legend.col <- c("grey", col[round(legend.y.col * length(col) / (maxbin.num - legend.min + 1)) - legend.min + 1])
		if(plot)	legend(x=(chorm.maxlen + chorm.maxlen/100), y=( -width/2.5 - band * (length(chr.num) - length(chr.num) - 1)), title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=legend.col,
			cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
		if(!plot)	return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
	}

    if (is(MVP, "list")) {
        MVP.res <- cbind(MVP$glm.results, MVP$mlm.results, MVP$farmcpu.results)
        Cnames <- colnames(MVP.res)
        Cnames <- Cnames[seq(3, ncol(MVP.res), 3)]
        MVP.res <- MVP.res[, seq(3, ncol(MVP.res), 3)]
        MVP.res[is.na(MVP.res)] <- 1
        MVP.res <- cbind(MVP$map[,c(1:3)], MVP.res)
        Pmap <- MVP.res
        colnames(Pmap)[-c(1:3)] <- Cnames
    } else {
        Pmap <- MVP
    }
    file = file.type
	if(sum(plot.type %in% "b")==1) plot.type=c("c","m","q","d")

	taxa=colnames(Pmap)[-c(1:3)]
	if(!is.null(memo) && memo != "")	memo <- paste("_", memo, sep="")
	if(length(taxa) == 0)	taxa <- paste("Col", 1:(ncol(Pmap)-3), sep="")
	taxa <- paste(taxa, memo, sep="")

	#SNP-Density plot
	if("d" %in% plot.type){
		if(verbose)	cat("SNP_Density Plotting...\n")
		if(file.output){
			ht=ifelse(is.null(height), 7, height)
			wh=ifelse(is.null(width), 9, width)
			if(file=="jpg")	jpeg(paste(outpath, "/SNP_Density.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
			if(file=="pdf")	pdf(paste(outpath, "/SNP_Density.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
			if(file=="tiff")	tiff(paste(outpath, "/SNP_Density.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
			par(xpd=TRUE)
		}else{
			ht=ifelse(is.null(height), 7, height)
			wh=ifelse(is.null(width), 9, width)
			if(is.null(dev.list()))	dev.new(width = wh,height=ht)
			par(xpd=TRUE)
		}
		if(!is.null(bin.range)){
			if(length(bin.range) != 2)	stop("Two values (min and max) should be provided for bin.range!")
			if(bin.range[1] == 0)	stop("Min value of bin.range should be more than 1!")
		}
		Densitplot(map=Pmap[,c(1:3)], chr.labels = chr.labels, col=chr.den.col, bin=bin.size, legend.min=bin.range[1], legend.max=bin.range[2], main=paste("The number of SNPs within ", bin.size/1e6, "Mb window size", sep=""))
		if(file.output)	dev.off()
	}

	if(length(plot.type) !=1 | (!"d" %in% plot.type)){
	
		#order Pmap by the name of SNP
		#Pmap=Pmap[order(Pmap[,1]),]
		Pmap <- as.matrix(Pmap)

		#scale and adjust the parameters
		cir.chr.h <- cir.chr.h/5
		cir.band <- cir.band/5
		if(!is.null(threshold)){
			threshold.col <- rep(threshold.col,length(threshold))
			threshold.lwd <- rep(threshold.lwd,length(threshold))
			threshold.lty <- rep(threshold.lty,length(threshold))
			signal.col <- rep(signal.col,length(threshold))
			signal.pch <- rep(signal.pch,length(threshold))
			signal.cex <- rep(signal.cex,length(threshold))
		}
		if(length(cex)!=3) cex <- rep(cex,3)
		if(!is.null(ylim)){
			if(length(ylim)==1) ylim <- c(0,ylim)
		}
		
		# if(is.null(conf.int.col))	conf.int.col <- NA
		# if(is.na(conf.int.col)){
		# 	conf.int=FALSE
		# }else{
		# 	conf.int=TRUE
		# }

		#get the number of traits
		R=ncol(Pmap)-3

		pch=rep(pch, R)

		#replace the non-euchromosome
		options(warn = -1)
		numeric.chr <- as.numeric(Pmap[, 2])
		options(warn = 0)
		max.chr <- max(numeric.chr, na.rm=TRUE)
		if(is.infinite(max.chr))	max.chr <- 0
		map.xy.index <- which(!numeric.chr %in% c(0:max.chr))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(Pmap[map.xy.index, 2])
			for(i in 1:length(chr.xy)){
				Pmap[Pmap[, 2] == chr.xy[i], 2] <- max.chr + i
			}
		}

		SNP_id <- Pmap[,1]

		#delete the column of SNPs names
		Pmap <- Pmap[,-1]

		Pmap <- matrix(as.numeric(Pmap), nrow(Pmap))
		order_index <- order(Pmap[, 1], Pmap[,2])

		#order the GWAS results by chromosome and position
		Pmap <- Pmap[order_index, ]
		SNP_id <- SNP_id[order_index]

		if(!is.null(highlight)){
			highlight <- as.character(as.matrix(highlight))
			highlight_index <- match(highlight, SNP_id)
			if(all(is.na(highlight_index)))	stop("No shared SNPs between Pmap and highlight!")
			highlight_index <- na.omit(highlight_index)
		}

		#get the index of chromosome
		chr <- unique(Pmap[,1])
		chr.ori <- chr
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
			}
		}

		pvalueT <- as.matrix(Pmap[,-c(1:2)])
		pvalue.pos <- Pmap[, 2]
		p0.index <- Pmap[, 1] == 0
		if(sum(p0.index) != 0){
			pvalue.pos[p0.index] <- 1:sum(p0.index)
		}
		pvalue.pos.list <- tapply(pvalue.pos, Pmap[, 1], list)
		
		#scale the space parameter between chromosomes
		if(!missing(band)){
			band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
		}else{
			band <- floor((sum(sapply(pvalue.pos.list, max))/100))
		}
		if(band==0)	band=1
		
		if(LOG10){
			pvalueT[pvalueT <= 0] <- NA
			pvalueT[pvalueT > 1] <- NA
		}
		Pmap[,-c(1:2)] <- pvalueT

		#set the colors for the plot
		#palette(heat.colors(1024)) #(heatmap)
		#T=floor(1024/max(pvalue))
		#plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
		if(is.vector(col)){
			col <- matrix(col,R,length(col),byrow=TRUE)
		}
		if(is.matrix(col)){
			#try to transform the colors into matrix for all traits
			col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=TRUE)
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
		
		#insert the space into chromosomes and return the midpoint of each chromosome
		ticks <- NULL
		pvalue.posN <- NULL
		#pvalue <- pvalueT[,j]
		if(Nchr == 1){
			bp <- ifelse((max_no_na(pvalue.pos.list[[1]]) -  min_no_na(pvalue.pos.list[[1]])) > 1000000, 1000000, 1000)
			bp_lab <- ifelse(bp == 1000000, " (Mb)", " (Kb)")
			pvalue.posN <- pvalue.pos.list[[1]] + band
			ticks <- seq(min_no_na(pvalue.pos.list[[1]]), max_no_na(pvalue.pos.list[[1]]), length=10)
			ticks <- seq(round(min_no_na(pvalue.pos.list[[1]]) / bp), round(max_no_na(pvalue.pos.list[[1]]) / bp), round((ticks[2]-ticks[1])/bp) + 0.5)
			if(!round(max_no_na(pvalue.pos.list[[1]]) / bp) %in% ticks){
				if(round(max_no_na(pvalue.pos.list[[1]]) / bp) - ticks[length(ticks)] > 0.5 * ticks[2])
				ticks <- c(ticks, round(max_no_na(pvalue.pos.list[[1]]) / bp))
			}
			ticks <- ticks[-1]
			chr.labels <- ticks
			ticks <- ticks * bp + band
		}else{
			for(i in 0:(Nchr-1)){
				if (i==0){
					#pvalue <- append(pvalue,rep(Inf,band),after=0)
					pvalue.posN <- pvalue.pos.list[[i+1]] + band
					ticks[i+1] <- max_no_na(pvalue.posN)-floor(max_no_na(pvalue.pos.list[[i+1]])/2)
				}else{
					#pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
					pvalue.posN <- c(pvalue.posN, max_no_na(pvalue.posN) + band + pvalue.pos.list[[i+1]])
					ticks[i+1] <- max_no_na(pvalue.posN)-floor(max_no_na(pvalue.pos.list[[i+1]])/2)
				}
			}
		}

		pvalue.posN.list <- tapply(pvalue.posN, Pmap[, 1], list)
		#NewP[[j]] <- pvalue
		
		#merge the pvalues of traits by column
		if(LOG10){
			logpvalueT <- -log10(pvalueT)
		}else{
			logpvalueT <- pvalueT
		}

		add <- list()
		for(i in 1:R){
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
		}

		TotalN <- max_no_na(pvalue.posN)

		if(length(chr.den.col) > 1){
			cir.density=TRUE
			den.fold <- 20
			density.list <- Densitplot(map=Pmap[,c(1,1,2)], col=chr.den.col, plot=FALSE, bin=bin.size, legend.min=bin.range[1], legend.max=bin.range[2])
			#list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
		}else{
			cir.density=FALSE
		}
		
		signal.line.index <- NULL
		if(!is.null(threshold)){
			if(!is.null(signal.line)){
				for(l in 1:R){
					if(LOG10){
						signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min_no_na(threshold)))
					}else{
						signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max_no_na(threshold)))
					}
				}
				signal.line.index <- unique(signal.line.index)
			}
		}
		signal.line.index <- pvalue.posN[signal.line.index]
	}
	
	#plot circle Manhattan
	if("c" %in% plot.type){
		#print("Starting Circular-Manhattan plot!",quote=F)
		if(file.output){
			ht=ifelse(is.null(height), 10, height)
			wh=ifelse(is.null(width), 10, width)
			if(file=="jpg")	jpeg(paste(outpath, "/Circular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
			if(file=="pdf")	pdf(paste(outpath, "/Circular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
			if(file=="tiff")	tiff(paste(outpath, "/Circular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
			par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
		}
		if(!file.output){
			ht=ifelse(is.null(height), 10, height)
			wh=ifelse(is.null(width), 10, width)
			if(is.null(dev.list()))	dev.new(width=wh, height=ht)
			par(pty="s", xpd=TRUE)
		}
		RR <- r+H*R+cir.band*R
		if(cir.density){
			plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
		}else{
			plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
		}
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
			
			if(verbose)	cat(paste("Circular_Manhattan Plotting ",taxa[i],"...\n",sep=""))
			pvalue <- pvalueT[,i]
			logpvalue <- logpvalueT[,i]
			if(is.null(ylim)){
				if(LOG10){
					Max <- max_ylim(-log10(min_no_na(pvalue)))
					Min <- 0
				}else{
					Max <- max_ylim(max_no_na(pvalue))
					#if(abs(Max)<=1)	Max <- max_no_na(pvalue)
					Min <- min_ylim(min_no_na(pvalue))
					#if(abs(Min)<=1)	Min <- min_no_na(pvalue)
				}
			}else{
				Max <- ylim[2]
				Min <- ylim[1]
			}
			Cpvalue <- (H*(logpvalue-Min))/(Max-Min)
			if(outward==TRUE){
				if(cir.chr==TRUE){
					
					#plot the boundary which represents the chromosomes
					polygon.num <- 1000
					for(k in 1:length(chr)){
						if(k==1){
							polygon.index <- seq(round(band/2)+1,-round(band/2)+max_no_na(pvalue.posN.list[[1]]), length=polygon.num)
							#change the axis from right angle into circle format
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
							}else{
								if(cir.density){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
								}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
								}
							}
						}else{
							polygon.index <- seq(1+round(band/2)+max_no_na(pvalue.posN.list[[k-1]]),-round(band/2)+max_no_na(pvalue.posN.list[[k]]), length=polygon.num)
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
							}else{
								if(cir.density){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
								}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
								}
							}		
						}
					}
					
					if(cir.density){

						segments(
							(RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							col=density.list$den.col, lwd=0.1
						)
						legend(
							x=RR+4*cir.chr.h,
							y=(RR+4*cir.chr.h)/2,
							title="", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
							cex=1, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=0.5, xjust=0, xpd=TRUE
						)
						
					}
					
					# XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					if(cir.density){
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
						circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
					}else{
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
						circle.plot(myr=RR,lwd=1.5,add=TRUE)
					}

				}
				
				X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
				Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
				points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				
				#plot the legend for each trait
				if(cir.legend==TRUE){
					#try to get the number after radix point
					if((Max-Min) > 1) {
						round.n=2
					}else{
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					text(-r/15,r+H*(i-0.94)+cir.band*(i-1),round(Min+(Max-Min)*0,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Min+(Max-Min)*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Min+(Max-Min)*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Min+(Max-Min)*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.06)+cir.band*(i-1),round(Min+(Max-Min)*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=ifelse(LOG10, H*(-log10(threshold[thr])-Min)/(Max-Min), H*(threshold[thr]-Min)/(Max-Min))
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								# warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
					}
				}
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify==TRUE){
							if(LOG10){
								threshold <- sort(threshold)
								significantline1=H*(-log10(max_no_na(threshold))-Min)/(Max-Min)
							}else{
								threshold <- sort(threshold, decreasing=TRUE)
								significantline1=H*(min_no_na(threshold)-Min)/(Max-Min)
							}
							
							p_amp.index <- which(Cpvalue>=significantline1)
							HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
										}else{
											significantline1=H*(threshold[ll]-Min)/(Max-Min)
										}
										p_amp.index <- which(Cpvalue>=significantline1)
										HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}else{
										if(LOG10){
											significantline0=H*(-log10(threshold[ll-1])-Min)/(Max-Min)
											significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
										}else{
											significantline0=H*(threshold[ll-1]-Min)/(Max-Min)
											significantline1=H*(threshold[ll]-Min)/(Max-Min)
										}
										p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
										HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}
								
									if(is.null(signal.col)){
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
									}else{
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[1],col=signal.col[ll])
									}
								}
						}
					}
				}

				if(!is.null(highlight)){
					points(X[highlight_index],Y[highlight_index],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
				}

				if(cir.chr==TRUE){
					ticks1=(RR+2*cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=(RR+2*cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						if(Nchr == 1){
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],paste(chr.labels[i], "Mb", sep=""),srt=angle,font=2,cex=cex.axis)
							}
						}else{
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
							}
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
						if(Nchr == 1){
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],paste(chr.labels[i], "Mb", sep=""),srt=angle,font=2,cex=cex.axis)
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
			if(outward==FALSE){
				if(cir.chr==TRUE){
					# XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)

					polygon.num <- 1000
					for(k in 1:length(chr)){
						if(k==1){
							polygon.index <- seq(round(band/2)+1,-round(band/2)+max_no_na(pvalue.posN.list[[1]]), length=polygon.num)
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
								if(is.null(chr.den.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									if(cir.density){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
									}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
									}
								}
						}else{
							polygon.index <- seq(1+round(band/2)+max_no_na(pvalue.posN.list[[k-1]]),-round(band/2)+max_no_na(pvalue.posN.list[[k]]), length=polygon.num)
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
							}else{
									if(cir.density){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
									}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
									}
							}	
						}
					}
					if(cir.density){

						segments(
							(RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							col=density.list$den.col, lwd=0.1
						)
						legend(
							x=RR+4*cir.chr.h,
							y=(RR+4*cir.chr.h)/2,
							title="", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
							cex=1, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=0.5, xjust=0, xpd=TRUE
						)
						
					}
					
					if(cir.density){
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
						circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
					}else{
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
						circle.plot(myr=RR,lwd=1.5,add=TRUE)
					}

				}
				
				X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
				Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
				points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				
				if(cir.legend==TRUE){
					
					#try to get the number after radix point
					if((Max-Min)<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					text(-r/15,r+H*(i-0.06)+cir.band*(i-1),round(Min+(Max-Min)*0,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Min+(Max-Min)*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Min+(Max-Min)*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Min+(Max-Min)*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.94)+cir.band*(i-1),round(Min+(Max-Min)*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
					
						for(thr in 1:length(threshold)){
							significantline1=ifelse(LOG10, H*(-log10(threshold[thr])-Min)/(Max-Min), H*(threshold[thr]-Min)/(Max-Min))
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								# warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						if(amplify==TRUE){
							if(LOG10){
								threshold <- sort(threshold)
								significantline1=H*(-log10(max_no_na(threshold))-Min)/(Max-Min)
							}else{
								threshold <- sort(threshold, decreasing=TRUE)
								significantline1=H*(min_no_na(threshold)-Min)/(Max-Min)
							}
							p_amp.index <- which(Cpvalue>=significantline1)
							HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
										}else{
											significantline1=H*(threshold[ll]-Min)/(Max-Min)
										}
										p_amp.index <- which(Cpvalue>=significantline1)
										HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}else{
										if(LOG10){
											significantline0=H*(-log10(threshold[ll-1])-Min)/(Max-Min)
											significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
										}else{
											significantline0=H*(threshold[ll-1]-Min)/(Max-Min)
											significantline1=H*(threshold[ll]-Min)/(Max-Min)
										}
										p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
										HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									
									}
								
									if(is.null(signal.col)){
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
									}else{
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[1],col=signal.col[ll])
									}
								}
						}
					}
				}
				
				if(!is.null(highlight)){
					points(X[highlight_index],Y[highlight_index],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
				}

				if(cir.chr==TRUE){
					ticks1=(RR+2*cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=(RR+2*cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						  angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						  text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						if(Nchr == 1){
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],paste(chr.labels[i], "Mb",sep=""),srt=angle,font=2,cex=cex.axis)
							}
						}else{
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
							}
						}
					}
				}else{
					ticks1=RR*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=RR*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						
							#adjust the angle of labels of circle plot
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						if(Nchr == 1){
							for(i in 1:length(ticks)){
								angle=360*(1-(ticks-round(band/2))[i]/TotalN)
								text(ticks1[i],ticks2[i],paste(chr.labels[i], "Mb",sep=""),srt=angle,font=2,cex=cex.axis)
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
		}
		if(file.output) dev.off()
		#print("Circular-Manhattan has been finished!",quote=F)
	}

	if("m" %in% plot.type){
		if(multracks==FALSE){
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			for(i in 1:R){
				colx=col[i,]
				colx=colx[!is.na(colx)]
				if(verbose)	cat(paste("Rectangular_Manhattan Plotting ",taxa[i],"...\n",sep=""))
					if(file.output){
						ht=ifelse(is.null(height), 6, height)
						wh=ifelse(is.null(width), 14, width)
						if(file=="jpg")	jpeg(paste(outpath, "/Rectangular-Manhattan.",taxa[i],".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
						if(file=="pdf")	pdf(paste(outpath, "/Rectangular-Manhattan.",taxa[i],".pdf",sep=""), width = wh,height=ht)
						if(file=="tiff")	tiff(paste(outpath, "/Rectangular-Manhattan.",taxa[i],".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
						par(mar = c(5,6,4,3),xaxs=xaxs,yaxs=yaxs,xpd=TRUE)
					}
					if(!file.output){
						ht=ifelse(is.null(height), 6, height)
						wh=ifelse(is.null(width), 14, width)
						if(is.null(dev.list()))	dev.new(width = wh, height = ht)
						par(xpd=TRUE)
					}
					
					pvalue=pvalueT[,i]
					logpvalue=logpvalueT[,i]
					if(is.null(ylim)){
						if(!is.null(threshold)){
							if(sum(threshold!=0)==length(threshold)){
								if(LOG10 == TRUE){
									Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),(-log10(min_no_na(threshold))))))
									Min <- 0
								}else{
									Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
									#if(abs(Max)<=1)	Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
									Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
									#if(abs(Min)<=1)	Min=min_no_na(c(min_no_na(pvalue),min_no_na(threshold)))
								}
							}else{
								if(LOG10){
									Max=max_ylim(-log10(min_no_na(pvalue)))
									Min<-0
								}else{
									Max=max_ylim(max_no_na(pvalue))
									#if(abs(Max)<=1)	Max=max_no_na(c(max_no_na(pvalue)))
									Min<-min_ylim(min_no_na(pvalue))
									#if(abs(Min)<=1)	Min=min_no_na(pvalue)
									# }else{
										# Max=max_no_na(ceiling(max_no_na(pvalue)))
									# }
								}
							}
						}else{
							if(LOG10){
									Max=max_ylim(-log10(min_no_na(pvalue)))
									Min<-0
							}else{
									Max=max_ylim(max_no_na(pvalue))
									#if(abs(Max)<=1)	Max=max_no_na(c(max_no_na(pvalue)))
									Min<-min_ylim(min_no_na(pvalue))
									#if(abs(Min)<=1)	Min=min_no_na(pvalue)
									# }else{
										# Max=max_no_na(ceiling(max_no_na(pvalue)))
									# }
							}
						}
						if((Max-Min)<=1){
							if(cir.density){
								plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ylab=ylab,
									cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
							}else{
								plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),ylab=ylab,
								cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
							}
						}else{
							if(cir.density){
								plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ylab=ylab,
								cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
							}else{
								plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),ylab=ylab,
								cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
							}
						}
					}else{
						Max <- max_no_na(ylim)
						Min <- min_no_na(ylim)
						if(cir.density){
							plot(pvalue.posN[logpvalue>=min_no_na(ylim)],logpvalue[logpvalue>=min_no_na(ylim)],pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim)],xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(min_no_na(ylim)-(Max-Min)/den.fold, max_no_na(ylim)),ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
						}else{
							plot(pvalue.posN[logpvalue>=min_no_na(ylim)],logpvalue[logpvalue>=min_no_na(ylim)],pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim)],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=ylim,ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
						}
					}
					# Max1 <- Max
					# Min1 <- Min
					# if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
					# if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
					if(!is.null(chr.labels)){
						if(Nchr == 1){
							axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
						}else{
							axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
							#axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
						}
					}else{
						axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
						#axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
					}
					axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
					if(is.null(ylim)){
						if((Max-Min)>1){
							axis(2, las=1, lwd=lwd.axis,cex.axis=cex.axis,font=2)
							axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
							legend.y <- tail(round(seq(Min,(Max),ceiling((Max-Min)/10)), 2), 1)
						}else{
							axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
							axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
							legend.y <- tail(seq(Min,Max,10^(-ceiling(-log10(max_no_na(abs(c(Max, Min))))))), 1)
						}
					}else{
						if(ylim[2]>1){
							axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
							axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
							legend.y <- tail(ylim[2], 1)
						}else{
							axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
							axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
							legend.y <- tail(ylim[2], 1)
						}
					}
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							for(thr in 1:length(threshold)){
								h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
								# print(h)
								# print(threshold.col[thr])
								# print(threshold.lty[thr])
								# print(threshold.lwd[thr])
								par(xpd=FALSE); abline(h=h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr]); par(xpd=TRUE)
							}
							if(amplify == TRUE){
								if(LOG10){
									threshold <- sort(threshold)
									sgline1=-log10(max_no_na(threshold))
								}else{
									threshold <- sort(threshold, decreasing=TRUE)
									sgline1=min_no_na(threshold)
								}

								sgindex=which(logpvalue>=sgline1)
								HY1=logpvalue[sgindex]
								HX1=pvalue.posN[sgindex]
								
								#cover the points that exceed the threshold with the color "white"
								points(HX1,HY1,pch=pch,cex=cex[2],col="white")
								
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											sgline1=-log10(threshold[ll])
										}else{
											sgline1=threshold[ll]
										}
										sgindex=which(logpvalue>=sgline1)
										HY1=logpvalue[sgindex]
										HX1=pvalue.posN[sgindex]
									}else{
										if(LOG10){
											sgline0=-log10(threshold[ll-1])
											sgline1=-log10(threshold[ll])
										}else{
											sgline0=threshold[ll-1]
											sgline1=threshold[ll]
										}
										sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
										HY1=logpvalue[sgindex]
										HX1=pvalue.posN[sgindex]
									}

									if(is.null(signal.col)){
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[2],col=rep(rep(colx,N[i]),add[[i]])[sgindex])
									}else{
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[2],col=signal.col[ll])
									}
									
								}
							}
						}
					}

					if(!is.null(highlight)){
						points(pvalue.posN[highlight_index],logpvalue[highlight_index],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
					}

					#if(!is.null(threshold) & !is.null(signal.line))	abline(v=pvalue.posN[which(pvalueT[,i] < min_no_na(threshold))],col="grey",lty=2,lwd=signal.line)
			
					if(is.null(ylim)){ymin <- Min}else{ymin <- min_no_na(ylim)}
					if(cir.density){
						for(yll in 1:length(pvalue.posN.list)){
							polygon(c(min_no_na(pvalue.posN.list[[yll]]), min_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]])), 
								c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
								ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
								col="grey", border="grey")
						}
						
						segments(
							pvalue.posN,
							ymin-0.5*(Max-Min)/den.fold,
							pvalue.posN,
							ymin-1.5*(Max-Min)/den.fold,
							col=density.list$den.col, lwd=0.1
						)
						legend(
							x=max_no_na(pvalue.posN)+band,
							y=legend.y,
							title="", legend=density.list$legend.y, pch=15, pt.cex = 2.5, col=density.list$legend.col,
							cex=0.8, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=1, xjust=0, xpd=TRUE
						)
						
					}
				if(box) box(lwd=lwd.axis)
				#if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))	abline(v=which(pvalueT[,i] < min_no_na(threshold)/max_no_na(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
				if(file.output)  dev.off()
			}
			#print("Rectangular-Manhattan has been finished!",quote=F)
		}else{
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			#print("Plotting in multiple tracks!",quote=F)
			if(file.output){
				ht=ifelse(is.null(height), 6, height)
				wh=ifelse(is.null(width), 14, width)
				if(file=="jpg")	jpeg(paste(outpath, "/Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi*R,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste(outpath, "/Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht*R)
				if(file=="tiff")	tiff(paste(outpath, "/Multracks.Rectangular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi*R,res=dpi)
				par(mfcol=c(R,1),mar=c(0+cex.lab, 6+cex.lab, 0, 2+cex.lab),oma=c(4,0,4,0),xaxs=xaxs,yaxs=yaxs,xpd=TRUE)
			}
			if(!file.output){
				ht=ifelse(is.null(height), 6, height)
				wh=ifelse(is.null(width), 14, width)
				if(is.null(dev.list()))	dev.new(width = wh, height = ht)
				par(xpd=TRUE)
			}
			for(i in 1:R){
				if(verbose)	cat(paste("Multracks_Rectangular Plotting ",taxa[i],"...\n",sep=""))
				colx=col[i,]
				colx=colx[!is.na(colx)]
				pvalue=pvalueT[,i]
				logpvalue=logpvalueT[,i]
				if(is.null(ylim)){
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							if(LOG10){
								Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(threshold)))))
								Min <- 0
							}else{
								Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
								#if(abs(Max)<=1)	Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
								Min<-min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
								#if(abs(Min)<=1)	Min=min_no_na(min_no_na(pvalue),min_no_na(threshold))
							}
						}else{
							if(LOG10){
								Max=max_ylim((-log10(min_no_na(pvalue))))
								Min<-0
							}else{
								Max=max_ylim((max_no_na(pvalue)))
								#if(abs(Max)<=1)	Max=max_no_na(max_no_na(pvalue))
								Min=min_ylim((min_no_na(pvalue)))
								#if(abs(Min)<=1)	Min=min_no_na(min_no_na(pvalue))
								# }else{
									# Max=max_no_na(ceiling(max_no_na(pvalue)))
								# }
							}	
						}
					}else{
						if(LOG10){
								Max=max_ylim((-log10(min_no_na(pvalue))))
								Min<-0
						}else{
								Max=max_ylim((max_no_na(pvalue)))
								#if(abs(Max)<=1)	Max=max_no_na(max_no_na(pvalue))
								Min=min_ylim((min_no_na(pvalue)))
								#if(abs(Min)<=1)	Min=min_no_na(min_no_na(pvalue))
								# }else{
									# Max=max_no_na(ceiling(max_no_na(pvalue)))
								# }
						}
					}
					if((Max-Min)<=1){
						plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=c(Min,Max),ylab=ylab,
							cex.axis=cex.axis*(R/2),cex.lab=cex.lab*(R/2),font=2,axes=FALSE,xlab="")
					}else{
						plot(pvalue.posN,logpvalue,pch=pch,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]]),xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=c(Min,Max),ylab=ylab,
							cex.axis=cex.axis*(R/2),cex.lab=cex.lab*(R/2),font=2,axes=FALSE,xlab="")
					}
				}else{
					Max <- max_no_na(ylim)
					Min <- min_no_na(ylim)
					plot(pvalue.posN[logpvalue>=min_no_na(ylim)],logpvalue[logpvalue>=min_no_na(ylim)],pch=pch,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim)],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=ylim,ylab=ylab,
						cex.axis=cex.axis*(R/2),cex.lab=cex.lab*(R/2),font=1,axes=FALSE,xlab="")
				}
				# Max1 <- Max
				# Min1 <- Min
				# if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
				# if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
				
				#add the names of traits on plot 
	
				text(max_no_na(pvalue.posN)*1.02,Max*0.5,labels=taxa[i],adj=0.5,font=1,cex=cex.lab*(R/2),srt=270)
		
				if(i == R){
					if(is.null(chr.labels)){
						axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.ori),padj=1)
						axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
					}else{
						if(Nchr == 1){
							axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
							axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
						}else{
							axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.labels))
							axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
						}
					}
				}
				#if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=2,cex=xn)
				if(is.null(ylim)){
					if((Max-Min)>1){
						axis(2, las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
						axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
					}else{
						axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
						axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
					}
				}else{
					if(ylim[2]>1){
						axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
						axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
					}else{
						axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
						axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
							par(xpd=FALSE); abline(h=h,col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr]); par(xpd=TRUE)
						}
						if(amplify==TRUE){
								if(LOG10){
									threshold <- sort(threshold)
									sgline1=-log10(max_no_na(threshold))
								}else{
									threshold <- sort(threshold, decreasing=TRUE)
									sgline1=min_no_na(threshold)
								}
								sgindex=which(logpvalue>=sgline1)
								HY1=logpvalue[sgindex]
								HX1=pvalue.posN[sgindex]
								
								#cover the points that exceed the threshold with the color "white"
								points(HX1,HY1,pch=pch,cex=cex[2]*R,col="white")
								
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											sgline1=-log10(threshold[ll])
										}else{
											sgline1=threshold[ll]
										}
										sgindex=which(logpvalue>=sgline1)
										HY1=logpvalue[sgindex]
										HX1=pvalue.posN[sgindex]
									}else{
										if(LOG10){
											sgline0=-log10(threshold[ll-1])
											sgline1=-log10(threshold[ll])
										}else{
											sgline0=threshold[ll-1]
											sgline1=threshold[ll]
										}
										sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
										HY1=logpvalue[sgindex]
										HX1=pvalue.posN[sgindex]
									}

									if(is.null(signal.col)){
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[2]*R,col=rep(rep(colx,N[i]),add[[i]])[sgindex])
									}else{
										points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*cex[2]*R,col=signal.col[ll])
									}
									
								}
						}
					}
				}

				if(!is.null(highlight)){
					points(pvalue.posN[highlight_index],logpvalue[highlight_index],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
				}

				#if(!is.null(threshold) & !is.null(signal.line))	abline(v=pvalue.posN[which(pvalueT[,i] < min_no_na(threshold))],col="grey",lty=2,lwd=signal.line)
			}
			
			#add the labels of X-axis
			#mtext(xlab,side=1,padj=2.5,font=2,cex=R*2/3)
			if(file.output) dev.off()
			
			if(file.output){
				ht=ifelse(is.null(height), 6, height)
				wh=ifelse(is.null(width), 14, width)
				if(file=="jpg")	jpeg(paste(outpath, "/Multraits.Rectangular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste(outpath, "/Multraits.Rectangular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
				if(file=="tiff")	tiff(paste(outpath, "/Multraits.Rectangular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
				par(mar = c(5,6,4,3),xaxs=xaxs,yaxs=yaxs,xpd=TRUE)
			}
			if(!file.output){
				ht=ifelse(is.null(height), 6, height)
				wh=ifelse(is.null(width), 14, width)
				if(is.null(dev.list()))	dev.new(width = wh, height = ht)
				par(xpd=TRUE)
			}
			
			pvalue <- as.vector(Pmap[,3:(R+2)])
			if(is.null(ylim)){
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(LOG10){
							Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(threshold)))))
							Min<-0
						}else{
							Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
							# if(abs(Max)<=1)	Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
							Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
							# if(abs(Min)<=1)	Min=min_no_na(c(min_no_na(pvalue),min_no_na(threshold)))
						}
					}else{
						if(LOG10){
							Max=max_ylim((-log10(min_no_na(pvalue))))
							Min=0
						}else{
							Max=max_ylim((max_no_na(pvalue)))
							# if(abs(Max)<=1)	Max=max_no_na(max_no_na(pvalue))
							Min<- min_ylim((min_no_na(pvalue)))
							# if(abs(Min)<=1)	Min=min_no_na(min_no_na(pvalue))
							# }else{
								# Max=max_no_na(ceiling(max_no_na(pvalue)))
							# }
						}	
					}
				}else{
					if(LOG10){
							Max=max_ylim((-log10(min_no_na(pvalue))))
							Min=0
					}else{
							Max=max_ylim((max_no_na(pvalue)))
							# if(abs(Max)<=1)	Max=max_no_na(max_no_na(pvalue))
							Min<- min_ylim((min_no_na(pvalue)))
							# if(abs(Min)<=1)	Min=min_no_na(min_no_na(pvalue))
							# }else{
								# Max=max_no_na(ceiling(max_no_na(pvalue)))
					}
				}
				if((Max-Min)<=1){
					if(cir.density){
						plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
					}else{
						plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
					}
				}else{
					if(cir.density){
						plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
					}else{
						plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),ylab=ylab,
							cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
					}
				}
			}else{
				Max <- max_no_na(ylim)
				Min <- min_no_na(ylim)
				if(cir.density){
					plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,1.01*max_no_na(pvalue.posN)),ylim=c(min_no_na(ylim)-Max/den.fold,Max),ylab=ylab,
						cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
				}else{
					plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=ylim,ylab=ylab,
						cex.axis=cex.axis,cex.lab=cex.lab,font=2,axes=FALSE,xlab=xlab,main="")
				}
			}
			# Max1 <- Max
			# Min1 <- Min
			# if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
			# if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
			if(!is.null(ylim)){
				legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,ylim[2]*1.2,taxa,col=t(col)[1:R],pch=pch,text.font=6,box.col=NA,horiz=TRUE,xjust=0.5)
			}else{
				legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max*1.2,taxa,col=t(col)[1:R],pch=pch,text.font=6,box.col=NA,horiz=TRUE,xjust=0.5)
			}
			if(is.null(chr.labels)){
				axis(1, at=c(min_no_na(pvalue.posN)-band,ticks),lwd=lwd.axis,cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
				axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
			}else{
				if(Nchr == 1){
					axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
					axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
				}else{
					axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
					axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
				}
			}
			if(is.null(ylim)){
				if((Max-Min)>1){
					#print(seq(0,(Max+1),ceiling((Max+1)/10)))
					axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
					axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
					legend.y <- tail(round(seq(0,(Max),ceiling((Max)/10)),2), 1)
				}else{
					axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
					axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
					legend.y <- tail(seq(0,Max,10^(-ceiling(-log10(max_no_na(abs(c(Max, Min))))))), 1)
				}
			}else{
				if(ylim[2]>1){
					axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
					axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
					legend.y <- tail(ylim[2], 1)
				}else{
					axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
					axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
					legend.y <- tail(ylim[2], 1)
				}
			}
			do <- TRUE
			sam.index <- list()
			for(l in 1:R){
				sam.index[[l]] <- 1:nrow(Pmap)
			}
			
			#change the sample number according to Pmap
			#sam.num <- ceiling(nrow(Pmap)/100)
			sam.num <- 1000
			cat_bar <- seq(1, 100, 1)
			while(do){
				for(i in 1:R){
					if(length(sam.index[[i]]) < sam.num){
						plot.index <- sam.index[[i]]
					}else{
						plot.index <- sample(sam.index[[i]], sam.num, replace=FALSE)
					}
					sam.index[[i]] <- sam.index[[i]][-which(sam.index[[i]] %in% plot.index)]
					logpvalue=logpvalueT[plot.index,i]
					if(!is.null(ylim)){indexx <- logpvalue>=min_no_na(ylim)}else{indexx <- 1:length(logpvalue)}
					points(pvalue.posN[plot.index][indexx],logpvalue[indexx],pch=pch[i],cex=cex[2],col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255))
					#if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))	abline(v=which(pvalueT[,i] < min_no_na(threshold)/max_no_na(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
				}
				if(verbose){
					progress <- round((nrow(Pmap) - length(sam.index[[i]])) * 100 / nrow(Pmap))
					if(progress %in% cat_bar){
						cat("Multraits_Rectangular Plotting...(finished ", progress, "%)\r", sep="")
						cat_bar <- cat_bar[cat_bar != progress]
						if(progress == 100)	cat("\n")
					}
				}
				if(length(sam.index[[i]]) == 0) do <- FALSE
			}
			
			# for(i in 1:R){
				# logpvalue=logpvalueT[,i]
				# points(pvalue.posN,logpvalue,pch=pch,cex=cex[2],col=t(col)[i])
			# }
			
			if(!is.null(threshold)){
				if(sum(threshold!=0)==length(threshold)){
					for(thr in 1:length(threshold)){
						h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
						par(xpd=FALSE); abline(h=h,col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr]); par(xpd=TRUE)
					}
				}
			}
			if(is.null(ylim)){ymin <- Min}else{ymin <- min_no_na(ylim)}
			if(cir.density){
						for(yll in 1:length(pvalue.posN.list)){
							polygon(c(min_no_na(pvalue.posN.list[[yll]]), min_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]])), 
								c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
								ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
								col="grey", border="grey")
						}
						
						segments(
							pvalue.posN,
							ymin-0.5*(Max-Min)/den.fold,
							pvalue.posN,
							ymin-1.5*(Max-Min)/den.fold,
							col=density.list$den.col, lwd=0.1
						)
						legend(
							x=max_no_na(pvalue.posN)+band,
							y=legend.y,
							title="", legend=density.list$legend.y, pch=15, pt.cex = 2.5, col=density.list$legend.col,
							cex=0.8, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=1, xjust=0, xpd=TRUE
						)
						
			}
			if(file.output) dev.off()
			
		}
	}
		
	if("q" %in% plot.type){
		#print("Starting QQ-plot!",quote=F)
		if(multracks){
			if(file.output){
				ht=ifelse(is.null(height), 5.5, height)
				wh=ifelse(is.null(width), 2.5, width)
				if(file=="jpg")	jpeg(paste(outpath, "/Multracks.QQplot.",paste(taxa,collapse="."),".jpg",sep=""), width = R*wh*dpi,height=ht*dpi,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste(outpath, "/Multracks.QQplot.",paste(taxa,collapse="."),".pdf",sep=""), width = R*wh,height=ht)
				if(file=="tiff")	tiff(paste(outpath, "/Multracks.QQplot.",paste(taxa,collapse="."),".tiff",sep=""), width = R*wh*dpi,height=ht*dpi,res=dpi)
				par(mfcol=c(1,R),mar = c(0,1,4,1.5),oma=c(3,5,0,0),xpd=TRUE)
			}else{
				ht=ifelse(is.null(height), 5.5, height)
				wh=ifelse(is.null(width), 2.5, width)
				if(is.null(dev.list()))	dev.new(width = wh*R, height = ht)
				par(xpd=TRUE)
			}
			log.Quantiles.max_no_na <- NULL
			for(i in 1:R){
				if(verbose)	cat(paste("Multracks_QQ Plotting ",taxa[i],"...\n",sep=""))		
				P.values=as.numeric(Pmap[,i+2])
				P.values=P.values[!is.na(P.values)]
				if(LOG10){
					P.values=P.values[P.values>0]
					P.values=P.values[P.values<1]
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				log.Quantiles.max_no_na <- c(log.Quantiles.max_no_na, max_no_na(log.Quantiles))
				if(LOG10){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}
				
				#calculate the confidence interval of QQ-plot
				if(conf.int){
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
				}else{
					c05 <- 1
					c95 <- 1
				}
				
				YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), floor(max_no_na(log.P.values)+1)))
				if(is.null(ylim)){
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,YlimMax),xlab ="", ylab="", main = taxa[i])
				}else{
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,max(ylim)),xlab ="", ylab="", main = taxa[i])
				}
				axis(1, at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				axis(2, las=1, lwd=lwd.axis,cex.axis=cex.axis)
				
				#plot the confidence interval of QQ-plot
				if(conf.int)	polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
				
				if(!is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1, col = threshold.col[1],lwd=2); par(xpd=TRUE)}
				points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						thre.line=-log10(min_no_na(threshold))
						if(amplify==TRUE){
							thre.index=which(log.P.values>=thre.line)
							if(length(thre.index)!=0){
							
								#cover the points that exceed the threshold with the color "white"
								points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
								if(is.null(signal.col)){
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch[1],cex=signal.cex[1])
								}else{
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
								}
							}
						}
					}
				}
			}
			if(box)	box(lwd=lwd.axis)
			if(file.output) dev.off()
			if(R > 1){
				signal.col <- NULL
				if(file.output){
					ht=ifelse(is.null(height), 5.5, height)
					wh=ifelse(is.null(width), 5.5, width)
					if(file=="jpg")	jpeg(paste(outpath, "/Multraits.QQplot.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
					if(file=="pdf")	pdf(paste(outpath, "/Multraits.QQplot.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
					if(file=="tiff")	tiff(paste(outpath, "/Multraits.QQplot.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
					par(mar = c(5,5,4,2),xpd=TRUE)
				}else{	
					ht=ifelse(is.null(height), 5.5, height)
					wh=ifelse(is.null(width), 5.5, width)
					dev.new(width = wh, height = ht)
					par(xpd=TRUE)
				}
				p_value_quantiles=(1:nrow(Pmap))/(nrow(Pmap)+1)
				log.Quantiles <- -log10(p_value_quantiles)
											
				# calculate the confidence interval of QQ-plot
				if(conf.int){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N1)
						if(xi==0)xi=1
						c95[j] <- qbeta(0.95,xi,N1-xi+1)
						c05[j] <- qbeta(0.05,xi,N1-xi+1)
					}
					index=length(c95):1
				}
				
				if(!conf.int){c05 <- 1; c95 <- 1}
				
				if(is.null(ylim)){
					Pmap.min_no_na <- Pmap[,3:(R+2)]
					YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), -log10(min_no_na(Pmap.min_no_na[Pmap.min_no_na > 0]))))
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0, floor(YlimMax+1)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = "QQplot")
				}else{
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0, max(ylim)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = "QQplot")
				}
				legend("topleft",taxa,col=t(col)[1:R],pch=19,text.font=6,box.col=NA)
				axis(1, at=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), cex.axis=cex.axis)
				axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis)
				
				for(i in 1:R){
					if(verbose)	cat(paste("Multraits_QQ Plotting ",taxa[i],"...\n",sep=""))
					P.values=as.numeric(Pmap[,i+2])
					P.values=P.values[!is.na(P.values)]
					if(LOG10){
						P.values=P.values[P.values>0]
						P.values=P.values[P.values<1]
						N=length(P.values)
						P.values=P.values[order(P.values)]
					}else{
						N=length(P.values)
						P.values=P.values[order(P.values,decreasing=TRUE)]
					}
					p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
					log.Quantiles <- -log10(p_value_quantiles)
					if(LOG10){
						log.P.values <- -log10(P.values)
					}else{
						log.P.values <- P.values
					}

					#calculate the confidence interval of QQ-plot
					if(conf.int){
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
					}else{
						c05 <- 1
						c95 <- 1
					}
	
					# plot the confidence interval of QQ-plot
					if(conf.int)	polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
				
						
					if((i == 1) & !is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1, col = threshold.col[1],lwd=2); par(xpd=TRUE)}
					points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
						
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							thre.line=-log10(min_no_na(threshold))
							if(amplify==TRUE){
								thre.index=which(log.P.values>=thre.line)
								if(length(thre.index)!=0){
								
									# cover the points that exceed the threshold with the color "white"
									points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
									if(is.null(signal.col)){
										points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch[1],cex=signal.cex[1])
									}else{
										points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
									}
								}
							}
						}
					}
				}
				if(box)	box(lwd=lwd.axis)
				if(file.output) dev.off()
			}
		}else{
			for(i in 1:R){
				if(verbose)	cat(paste("Q_Q Plotting ",taxa[i],"...\n",sep=""))
				if(file.output){
					ht=ifelse(is.null(height), 5.5, height)
					wh=ifelse(is.null(width), 5.5, width)
					if(file=="jpg")	jpeg(paste(outpath, "/QQplot.",taxa[i],".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
					if(file=="pdf")	pdf(paste(outpath, "/QQplot.",taxa[i],".pdf",sep=""), width = wh,height=ht)
					if(file=="tiff")	tiff(paste(outpath, "/QQplot.",taxa[i],".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
					par(mar = c(5,5,4,2),xpd=TRUE)
				}else{
					ht=ifelse(is.null(height), 5.5, height)
					wh=ifelse(is.null(width), 5.5, width)
					if(is.null(dev.list()))	dev.new(width = wh, height = ht)
					par(xpd=TRUE)
				}
				P.values=as.numeric(Pmap[,i+2])
				P.values=P.values[!is.na(P.values)]
				if(LOG10){
					P.values=P.values[P.values>0]
					P.values=P.values[P.values<1]
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				if(LOG10){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}
				
				#calculate the confidence interval of QQ-plot
				if(conf.int){
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
				}else{
					c05 <- 1
					c95 <- 1
				}
				if(is.null(ylim)){
					YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), floor(max_no_na(log.P.values)+1)))
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,YlimMax),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))
				}else{
					plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,max(ylim)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))		
				}
				axis(1, at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis)
				
				#plot the confidence interval of QQ-plot
				if(conf.int)	polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
				
				if(!is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1, col = threshold.col[1],lwd=2); par(xpd=TRUE)}
				points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						thre.line=-log10(min_no_na(threshold))
						if(amplify==TRUE){
							thre.index=which(log.P.values>=thre.line)
							if(length(thre.index)!=0){
							
								#cover the points that exceed the threshold with the color "white"
								points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
								if(is.null(signal.col)){
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch[1],cex=signal.cex[1])
								}else{
									points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
								}
							}
						}
					}
				}
				if(box)	box(lwd=lwd.axis)
				if(file.output) dev.off()
			}
		}
	}
	if(file.output & verbose)	cat(paste("Plots are stored in: ", getwd(), sep=""), "\n")
}

Densitplot <- function(map, col = c("darkgreen", "yellow", "red"), main = "SNP Density", bin = 1e6,
                       band = 3, width = 5, legend.len = 10, legend.max = NULL, legend.pt.cex = 3,
                       legend.cex = 1, legend.x.intersp = 1, legend.y.intersp = 1, plot = TRUE) {
    
    ## Step 1: preprocess map
    # filter map
    map <- as.matrix(map)
    map <- map[!is.na(map[, 2]), ]
    map <- map[!is.na(map[, 3]), ]
    map <- map[map[, 2] != 0, ]

    # get the number of chromosomes
    max.chr <- max(as.numeric(map[, 2]), na.rm = TRUE)
    if (is.infinite(max.chr)) { max.chr <- 0 }
    
    # deal with x,y
    map.xy.index <- suppressWarnings(which(!as.numeric(map[, 2]) %in% c(0:max.chr)))
    if (length(map.xy.index) != 0) {
        chr.xy <- unique(map[map.xy.index, 2])
        for (i in 1:length(chr.xy)) {
            map[map[, 2] == chr.xy[i], 2] <- max.chr + i
        }
    }
    
    # sort map
    map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
    chr <- as.numeric(map[, 2])
    pos <- as.numeric(map[, 3])
    chr.num <- unique(chr)
    chorm.maxlen <- max(pos)
    
    # Step 2: count SNP
    pos.x      <- list()
    col.index  <- list()
    maxbin.num <- NULL
    for (i in 1:length(chr.num)) {
        pos.x[[i]] <- pos[which(chr == chr.num[i])]
        cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
        if (cut.len <= 1) {
            col.index[[i]] = 1
        } else {
            cut.r          <- cut(pos.x[[i]], cut.len, labels = FALSE)
            eachbin.num    <- table(cut.r)
            maxbin.num     <- c(maxbin.num, max(eachbin.num))
            col.index[[i]] <- rep(eachbin.num, eachbin.num)
        }
    }
    Maxbin.num <- max(maxbin.num)
    maxbin.num <- Maxbin.num
    if (!is.null(legend.max)) {
        maxbin.num <- legend.max
    }
    col = colorRampPalette(col)(maxbin.num)
    col.seg = NULL
    for (i in 1:length(chr.num)) {
        if (!is.null(legend.max)) {
            if (legend.max < Maxbin.num) {
                col.index[[i]][col.index[[i]] > legend.max] <- legend.max
            }
        }
        col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
    }
    if (length(map.xy.index) != 0) {
        for (i in 1:length(chr.xy)) {
            chr.num[chr.num == max.chr + i] <- chr.xy[i]
        }
    }
    chr.num <- rev(chr.num)

    # image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
    # round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
    # t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
    
    ## Step 3: Deal with legend label and color.
    len      <- round(seq(0, maxbin.num, length = legend.len))[2]
    legend.y <- seq(0, maxbin.num, len)

    if (!is.null(legend.max)) {
        if (legend.max < Maxbin.num) {
            if (!maxbin.num %in% legend.y) {
                legend.y <- c(legend.y, paste0(">=", maxbin.num))
            } else {
                legend.y[length(legend.y)] <- paste0(">=", maxbin.num)
            }
            legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
        } else {
            if (!maxbin.num %in% legend.y) { 
                legend.y <- c(legend.y, maxbin.num) 
            }
            legend.y.col <- c(legend.y[-1])
        }
    } else {
        if (!maxbin.num %in% legend.y) {
            legend.y     <- c(legend.y, paste0(">", max(legend.y)))
            legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
        } else {
            legend.y.col <- c(legend.y[-1])
        }
    }
    
    legend.y.col <- as.numeric(legend.y.col)
    legend.col   <- c("grey", col[round(legend.y.col * length(col) / maxbin.num)])
    
    ## Step 4: Plot or Return
    if (plot) {
        # draw a canvas
        plot(
            NULL,
            xlim = c(0, chorm.maxlen + chorm.maxlen / 10),
            ylim = c(0, length(chr.num) * band + band),
            main = main,
            axes = FALSE,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i"
        )
        
        # draw each Chr
        for (i in 1:length(chr.num)) {
            # draw bg
            polygon(
                x = c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])),
                y = c(
                    -width / 5 - band * (i - length(chr.num) - 1),
                     width / 5 - band * (i - length(chr.num) - 1),
                     width / 5 - band * (i - length(chr.num) - 1),
                    -width / 5 - band * (i - length(chr.num) - 1)
                ),
                col = "grey",
                border = "grey"
            )
            # draw fg
            segments(
                x0 = pos.x[[i]],
                y0 = -width / 5 - band * (i - length(chr.num) - 1),
                x1 = pos.x[[i]],
                y1 = width / 5 - band * (i - length(chr.num) - 1),
                col = col[round(col.index[[i]] * length(col) / maxbin.num)],
                lwd = 1
            )
        }
        # draw Chr label
        mtext(
            at = seq(band, length(chr.num) * band, band),
            text = paste0("Chr", chr.num),
            side = 2,
            las = 2,
            font = 1,
            cex = 0.6,
            line = 0.2
        )
        
        # draw physical distance ruler
        axis(
            side = 3,
            at = seq(0, chorm.maxlen, length = 10),
            labels = c(NA, paste0(round((
                seq(0, chorm.maxlen, length = 10)
            )[-1] / 1e6, 0), "Mb")),
            font = 1,
            cex.axis = 0.8,
            tck = 0.01,
            lwd = 2,
            padj = 1.2
        )
        
        # draw legend
        legend(
            x = (chorm.maxlen + chorm.maxlen / 100),
            y = (-width / 2.5 - band * (length(chr.num) - length(chr.num) - 1)),
            legend = legend.y,
            col = legend.col,
            pch = 15,
            bty = "n",
            cex = legend.cex,
            pt.cex = legend.pt.cex,
            xjust = 0,
            yjust = 0,
            x.intersp = legend.x.intersp,
            y.intersp = legend.y.intersp,
            title = "",
            xpd = TRUE
        )
    } else {
        return(list(den.col = col.seg, legend.col = legend.col, legend.y = legend.y))
    }
}

#' SNP Density
#'
#' @param Pmap P value Map
#' @param col The color vector
#' @param dpi Number. Dots per inch for .jpg and .tiff files
#' @param bin.size the window size for counting SNP number
#' @param bin.max maximum SNP number, for winows, which has more SNPs than bin.max, will be painted in same color
#' @param memo Character. A text marker on output files
#' @param outpath Only when file.output = TRUE, determines the path of the output file
#' @param file.type format of output figure
#' @param file.output Whether to output the file
#' @param verbose whether to print detail.
#' 
#' @export
#' @return 
#' Output file:
#' <memo>.SNP_Density.<type>
#'
#' @examples
#' data(pig60K, package = "rMVP")
#' 
#' MVP.Report.Density(pig60K, file.output=FALSE)
#' 
MVP.Report.Density <- function(Pmap, col = c("darkgreen", "yellow", "red"), dpi = 300, outpath = getwd(), memo = 'MVP',
                               bin.size = 1e6, bin.max = NULL, file.type = "jpg", file.output = TRUE, verbose = TRUE) {
    logging.log("SNP_Density Plotting", "\n", verbose = verbose)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    
    w <- 9
    h <- 7
    if (file.output) {
        name <- file.path(outpath, paste0(memo, ".SNP_Density"))
        switch(file.type,
               jpg = jpeg(paste0(name, ".jpg"), width = w * dpi, height = h * dpi, res = dpi, quality = 100),
               pdf = pdf(paste0(name, ".pdf"), width = w, height = h),
               tiff = tiff(paste0(name, ".tiff"), width = w * dpi, height = h * dpi, res = dpi)
        )
        par(xpd = TRUE)
    } else {
        if (is.null(dev.list())) { dev.new(width = w, height = h) }
        par(xpd = TRUE)
    }
    
        Densitplot(map = Pmap[,c(1:3)], col = col, bin = bin.size, legend.max = bin.max, main = paste0("The number of SNPs within ", bin.size / 1e6, "Mb window size"))
    
    if (file.output) { dev.off() }
}

filter.points <- function(x, y, w, h, dpi=300, scale=1) {
    x <- ceiling((x - min(x)) / (max(x) - min(x)) * w * dpi / scale)
    y <- ceiling((y - min(y)) / (max(y) - min(y)) * h * dpi / scale)
    index <- !duplicated(cbind(x, y))
}

#' QQ Plot
#'
#' @param P.values P values
#' @param taxa_name The identifier of the phenotype will be used to generate a
#'    portion of the image file name. If the title parameter is NULL, it will 
#'    also be part of the title.
#' @param col default color is "blue"
#' @param cex A numerical value giving the amount by which plotting text and
#'    symbols should be magnified relative to the default. This starts as 1 
#'    when a device is opened, and is reset when the layout is changed, e.g. 
#'    by setting mfrow. see \code{\link[graphics]{par}}. 
#' @param threshold Number or Vector. The cutoff line on Manhattan plot, 
#'    e.g. Bonfferoni correction. More than one significant line can be added
#'    onto one figure. If threshold=0 or NULL, the threshold line will not be
#'    added.
#' @param amplify Logical value. If TRUE, the points that passed the threshold line will be highlighted
#' @param signal.cex Number. If "amplify" is TRUE, "signal.cex" is used to set the size of significant points
#' @param signal.pch Number. If "amplify" is TRUE, users can set the type of significant points
#' @param signal.col Character. If "amplify" is TRUE, "signal.col" is used to set the color of significant points, if "signal.col" is NULL, the colors of significant points will not be changed
#' @param conf.int Whether to draw a confidence interval
#' @param cex.axis a number, controls the size of numbers of X-axis and the size of labels of circle plot.
#' @param conf.int.col a character, the color of the confidence interval on QQ-plot.
#' @param threshold.col Character or Vector. The colors of threshold lines
#' @param outpath Only when file.output = TRUE, determines the path of the output file
#' @param file.type A string or NULL is used to determine the type of output 
#'    file. Can be "jpg", "pdf", "tiff". If it is NULL, it will use 
#'    \code{\link[grDevices:dev]{dev.new()}} to create a new graphics device 
#'    in the current environment, which may be RStudioGD or the default 
#'    device of the system.
#' @param memo the prefix of the output image file.
#' @param box A Boolean value that controls whether to draw a box around
#'    QQplot.
#' @param dpi a number, the picture element for .jpg and .tiff files. The default is 300.
#' @param verbose whether to print detail.
#' @param file.output Logical value. If TRUE, the figures will be generated. 
#' 
#' @return
#' Output file:
#' <memo>.QQplot.<taxa_name>.<type>
#' @export
#'
#' @examples
#' data(pig60K, package = "rMVP")
#' 
#' MVP.Report(pig60K, plot.type="q", file.output=FALSE)
#' 
MVP.Report.QQplot <-
    function(P.values,
             taxa_name,
             col=c("blue"),
             cex=0.5,
             threshold=NULL,
             amplify=TRUE,
             signal.col="red",
             signal.pch=19,
             signal.cex=0.8,
             conf.int=TRUE,
             cex.axis=1,
             conf.int.col="grey",
             threshold.col="red",
             outpath=getwd(),
             file.type="jpg",
             memo="MVP",
             box=TRUE,
             dpi=300,
             file.output=TRUE,
             verbose=TRUE
) {
    logging.log(paste0("Q_Q Plotting ", taxa_name), "\n", verbose = verbose)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    w <- 5.5
    h <- 5.5
    
    # setup output device
    if (file.output && !is.null(file.type)) {
        name <- file.path(outpath, paste0(memo, ".QQplot.", taxa_name))
        switch(file.type,
            jpg = jpeg(paste0(name, ".jpg"), width = w * dpi, height = h * dpi, res = dpi, quality = 100),
            pdf = pdf(paste0(name, ".pdf"), width = w, height = h),
            tiff = tiff(paste0(name, ".tiff"), width = w * dpi, height = h * dpi, res = dpi)
        )
        par(mar = c(5, 5, 4, 2), xpd = TRUE)
    } else {
        if (is.null(dev.list())) { dev.new(width = w, height = h) }
        par(xpd = TRUE)
    }
    
    # filter & order P-values
    P.values <- as.numeric(P.values)
    P.values <- P.values[!is.na(P.values)]
    P.values <- P.values[P.values > 0]
    P.values <- P.values[P.values <= 1]
    P.values <- -log10(P.values[order(P.values)])
    
    N <- length(P.values)
    quantiles <- -log10((1:N) / (N + 1))
    
    # filter
    is_visable <- filter.points(quantiles, P.values, w, h, dpi)
    
    # calculate the confidence interval of QQ-plot
    c95 <- 1
    c05 <- 1
    if (conf.int) {
        vi <- which(is_visable)
        c95 <- rep(NA, length(vi))
        c05 <- rep(NA, length(vi))
        for (i in 1:length(vi)) {
            xi <- ceiling((10 ^ -quantiles[vi[i]]) * N)
            if (xi == 0) { xi <- 1 }
            c95[i] <- qbeta(0.95, xi, N - xi + 1)
            c05[i] <- qbeta(0.05, xi, N - xi + 1)
        }
        c95 <- -log10(c95)
        c05 <- -log10(c05)
    }
    XlimMax <- ceiling(max(quantiles))
    YlimMax <- ceiling(max(c05, c95, P.values))
    
    plot(
        NULL,
        xlim = c(0, XlimMax),
        ylim = c(0, YlimMax),
        axes = FALSE,
        cex.axis = cex.axis,
        cex.lab = 1.2,
        xlab = expression(Expected~~-log[10](italic(p))),
        ylab = expression(Observed~~-log[10](italic(p))),
        main = paste("QQplot of", taxa_name)
    )
    axis(
        1,
        at = seq(0, XlimMax, ceiling(XlimMax / 10)),
        labels = seq(0, XlimMax, ceiling(XlimMax / 10)),
        cex.axis = cex.axis
    )
    axis(
        2,
        at = seq(0, YlimMax, ceiling(YlimMax / 10)),
        labels = seq(0, YlimMax, ceiling(YlimMax / 10)),
        cex.axis = cex.axis
    )
    
    # plot the confidence interval of QQ-plot
    if (conf.int) {
        polygon(
            x = c(rev(quantiles[is_visable]), quantiles[is_visable]),
            y = c(rev(c05), c95),
            col = conf.int.col,
            border = conf.int.col
        )
    }
    
    if (!is.null(threshold.col)) {
        abline(a = 0, b = 1, col = threshold.col, lwd = 2, xpd = FALSE)
    }
    
    if (is.null(signal.col)) { signal.col <- col }
    is_signal <- rep(FALSE, N)
    if (!is.null(threshold) && !(0 %in% threshold) && amplify) {
        is_signal <- (P.values >= -log10(min(threshold)))
    }
    
    points(quantiles[is_visable && !is_signal], P.values[is_visable && !is_signal], col = col, pch = 19, cex = cex)
    points(quantiles[is_visable && is_signal], P.values[is_visable && is_signal], col = signal.col, pch = signal.pch, cex = signal.cex)
    
    if (box) { box() }
    if (file.output && !is.null(file.type)) { dev.off() }
}

#' Phenotype distribution histogram
#'
#' @param phe phenotype data
#' @param col The color vector of the histogram. If the number of colors is 
#'    less than break.n, the color will be reused. If the number of colors is 
#'    greater than break.n, only the previous break.n colors will be used.
#' @param breakNum the number of cells for the histogram. The default value 
#'    is 15.
#' @param test.method The method used to test the normal distribution. The options are "auto", 
#'    "Shapiro-Wilk", "Kolmogorov-Smirnov", and NULL. When set to "auto", "Shapiro- Wilk" method,
#'    "Kolmogorov-Smirnov" method will be used when it is greater than 5000, and it will not be 
#'    tested when set to NULL.
#' @param outpath Effective only when file.output = TRUE, determines the path of the output file
#' @param file.type A string or NULL is used to determine the type of output 
#'    file. Can be "jpg", "pdf", "tiff". If it is NULL, it will use 
#'    \code{\link[grDevices:dev]{dev.new()}} to create a new graphics device 
#'    in the current environment, which may be RStudioGD or the default 
#'    device of the system.
#' @param dpi The resolution of the image, specifying how many pixels 
#'    per inch.
#' @param file.output Logical value. If TRUE, the figures will be generated. 
#' @param memo Character. A text marker on output files
#' 
#' @export
#' @return
#' Output file:
#' MVP.Phe_Distribution.<trait>.<type>
#'
#' @examples
#' phePath <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phe <- read.table(phePath, header=TRUE)
#' 
#' MVP.Hist(phe, file.output = FALSE)
#' 
MVP.Hist <-
    function(phe,
             col = c("dodgerblue4",
                     "olivedrab4",
                     "violetred",
                     "darkgoldenrod1",
                     "purple4"),
             breakNum = 15,
             memo = TRUE,
             outpath = getwd(), 
             test.method = "auto",
             file.type = "pdf",
             file.output = TRUE,
             dpi = 300) {
    
    w <- 6
    h <- 6
    phex <- phe
    
    for (i in 2:ncol(phe)) {
        # create file
        trait <- colnames(phe)[i]
        name  <- file.path(outpath, paste0(memo, ".Phe_Distribution.", paste(trait, collapse = "_")))
        
        if (file.output && !is.null(file.type)) {
            if (file.type == "jpg") { jpeg(paste0(name, ".jpg"), width = w * dpi, height = h * dpi, res = dpi, quality = 100) }
            if (file.type == "pdf") { pdf(paste0(name, ".pdf"), width = w, height = h) }
            if (file.type == "tiff") { tiff(paste0(name, ".tiff"), width = w * dpi, height = h * dpi, res = dpi) }
        } else {
            if (is.null(dev.list())) { dev.new(width = w, height = h) }
        }
        
        
        phe    <- phe[!is.na(phe[, i]), ]
        Breaks <- seq(min(phe[, i], na.rm = TRUE), max(phe[, i], na.rm = TRUE), length = breakNum)
        
        # hist
        xx <-
            suppressWarnings(hist(
                x = phe[, i],
                plot = FALSE,
                breaks = Breaks,
                xlab = "",
                ylab = "Density",
                freq = FALSE,
                col = colorRampPalette(col)(breakNum),
                font = 2,
                font.lab = 2,
                main = paste0("Distribution of ", trait)
            ))
        
        maxY <- max(max(xx$density),  max(density(phe[, i])$y))
        
        hist(
            x = phe[, i],
            breaks = Breaks,
            xlab = "",
            ylab = "Density",
            ylim = c(0, maxY),
            freq = FALSE,
            col = colorRampPalette(col)(breakNum),
            font = 2,
            font.lab = 2,
            main = paste0("Distribution of ", trait)
        )
        
        lines(density(phe[, i]), lwd = 2)
        
        # test normal distribution
        if (is.null(test.method)) {
            if (test.method == "auto") {
                if (length(phe[, i]) <= 5000) {
                    test.method <- "Shapiro-Wilk"
                } else {
                    test.method <- "Kolmogorov-Smirnov"
                }
            }
            
            if (test.method == "Shapiro-Wilk") {
                norm.p <- round(shapiro.test(phe[, i])$p, 4)
            } else if (test.method == "Kolmogorov-Smirnov") {
                norm.p <- round(ks.test(phe[, i], "pnorm")$p, 4)
            }
            
            text(xx$breaks[1], y = maxY*0.85, labels = paste0(test.method, ": ", norm.p), font = 2, adj = 0)
        }
        
        # draw text
        text(xx$breaks[1], y = maxY*0.95, labels = paste0("Mean: ", round(mean(phe[, i]), 2)), font = 2, adj = 0)
        text(xx$breaks[1], y = maxY*0.90, labels = paste0("  Sd: ", round(sd(phe[, i]), 2)), font = 2, adj = 0)
        if (file.output && !is.null(file.type)) { dev.off() }
        phe <- phex
    }
}

#' PCA Plot
#'
#' @param PCA Principal component analysis result, 2-column matrix
#' @param col colors for each cluster
#' @param pch Either an integer specifying a symbol or a single character to be 
#'    used as the default in plotting points. See \code{\link[graphics]{points}} 
#'    for possible values and their interpretation. Note that only integers and
#'    single-character strings can be set as a graphics parameter (and not NA 
#'    nor NULL).
#' @param class the class of all individuals, for example: "breed", "location"
#' @param legend.pos position of legend. default is "topright"
#' @param Ncluster cluster number
#' @param plot3D (DEPRECATED)if TRUE, plot PC figure in 3D format, it can be only used in windows and mac operation system, "rgl" package should be installed beforehead
#' @param file.type Character. Options are jpg, pdf, and tiff
#' @param dpi Number. Dots per inch for .jpg and .tiff files
#' @param box Logical value. If TRUE, the border line of Manhattan plot will be added
#' @param memo the prefix of the output image file.
#' @param file.output Logical value. If TRUE, the figures will be generated. 
#' @param outpath Effective only when file.output = TRUE, determines the path of the output file
#' @param verbose whether to print detail.
#' 
#' @export
#' @return
#' Output file:
#' MVP.PCA_2D.<type>
#'
#' @examples
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' geno <- attach.big.matrix(genoPath)
#' pca <- MVP.PCA(M=geno)
#' 
#' MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, 
#'     col=c("red", "green", "yellow"), file.output=FALSE, pch=19)
#' 
MVP.PCAplot <- function(PCA,
                        memo = "MVP",
                        col = NULL,
                        pch = NULL,
                        class = NULL,
                        legend.pos = "topright",
                        Ncluster = 1,
                        plot3D = FALSE,
                        file.type = "pdf",
                        dpi = 300,
                        box = FALSE,
                        file.output = TRUE,
                        outpath = getwd(),
                        verbose = TRUE
) {
    w <- 5.5
    h <- 5.5
    
    if(!is.null(class)){if(length(class) != nrow(PCA)) stop("the length of 'class' differs from the row of 'PCA'"); Ncluster <- length(unique(class))}
    if(!is.null(col)){
        col <- rep(col, Ncluster)
    }else{
        col <- c("dodgerblue1", "olivedrab3", "darkgoldenrod1", "red")
        col <- colorRampPalette(col)(Ncluster)
    }
    if(!is.null(pch)){
        pch <- rep(pch, Ncluster)
    }else{
        pch=1:Ncluster
    } 
    logging.log("PCA plot2d", "\n", verbose = verbose)
    
    if (file.output && !is.null(file.type)) {
        name <- file.path(outpath, paste0(memo, ".PCA_2D"))
        switch(file.type,
               jpg = jpeg(paste0(name, ".jpg"), width = w * dpi, height = h * dpi, res = dpi, quality = 100),
               pdf = pdf(paste0(name, ".pdf"), width = w, height = h),
               tiff = tiff(paste0(name, ".tiff"), width = w * dpi, height = h * dpi, res = dpi)
        )
    } else {
        if (is.null(dev.list())) { dev.new(width = w, height = h) }
    }

    if(is.null(class)){
        kc <- kmeans(PCA[,c(1,2)], Ncluster)
    }else{
        kc <- list(cluster=as.numeric(as.factor(class)))
    }
    plot(PCA[,1],PCA[,2],pch=pch[kc$cluster],col=col[kc$cluster],font=2,font.lab=2,xlab="PC1",ylab="PC2",axes=FALSE)
    axis(1, font.axis=2)
    if(!is.null(class)){
        legend(legend.pos,levels(as.factor(class)),col=col[1:Ncluster],pch=pch[1:Ncluster],text.font=6,box.col=NA)
    }
    axis(2, font.axis=2)
    if(box) box()
    if (file.output && !is.null(file.type)) { dev.off() }
    if(plot3D){
        logging.log("The 3D PCA map has been temporarily disabled, and we will improve this feature in subsequent versions.", "\n", verbose = verbose)
        # print("PCA plot3d")
        # par3d(cex=0.8,windowRect=c(100,100,600,600),font=4,userMatrix=matrix(c(0.88,0.47,-0.07,0,-0.14,0.40,0.90,0,0.45,-0.78,0.42,0,0,0,0,1),4,byrow=TRUE))
        # if(is.null(class)) kc <- kmeans(PCA, Ncluster)
        # plot3d(PCA,col=col[kc$cluster],type="s",size=1.5,top=FALSE)
        # #rgl.postscript("MVP.PCA_3D.pdf","pdf")
        # dir.create('PCA plot3D')
        # setwd('./PCA plot3D')
        # M <- par3d("userMatrix")
        # try(movie3d(par3dinterp(time=(0:2)*0.75,userMatrix=list(M,rotate3d(M, pi, 1, 0, 0), rotate3d(M, pi, 0, 1, 0))),fps=30,duration=3,movie="MVP.3D",dir=getwd(),top=FALSE), silent=TRUE)
        # movie3d(spin3d(rpm=20),duration=3,fps=15,frames="MVP.3D",convert=FALSE,dir=getwd(),top=FALSE)
        # setwd("../")
    }
}

