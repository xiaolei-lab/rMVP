MVP.Version <-
function(start=TRUE){
##############################################################################################
# MVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS
# Designed by Xiaolei Liu
# Writen by Lilin Yin
# Build date: Aug 30, 2017
# Last update: May 25, 2017
##############################################################################################
version <- "1.0.1"
authors <- "Lilin Yin, Haohao Zhang, Zhiwu Zhang, Xinyun Li, Xiaohui Yuan, Shuhong Zhao, and Xiaolei Liu"
title <- "A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS"
mail <- "xiaoleiliu@mail.hzau.edu.cn"
Total.len <- nchar(authors)+11
if(start){
cat(paste("#", paste(rep("-", floor((Total.len-14)/2)), collapse=""), "Welcome to MVP", paste(rep("-", ceiling((Total.len-14)/2)), collapse=""), "#",sep=""), "\n")
cat("#", paste(rep(" ", floor((Total.len-nchar(title))/2)-2), collapse=""), title, paste(rep(" ", ceiling((Total.len-nchar(title))/2-2)), collapse=""), "#", "\n")
cat(paste("# ", "Version: ", version, paste(rep(" ", Total.len-nchar(version)-10), collapse=""), "#", sep=""), "\n")
cat(paste("#", "Authors:", authors,"#", "\n"))
cat(paste("#", "Contact:", mail, paste(rep(" ", Total.len-nchar(mail)-12), collapse=""), "#", "\n"))
cat(paste("#",  paste(rep("-", Total.len), collapse=""), "#", sep=""), "\n")
}else{
cat(paste("#",  paste(rep("-", floor((Total.len-16)/2)), collapse=""), "MVP ACCOMPLISHED", paste(rep("-", ceiling((Total.len-16)/2)), collapse=""), "#", sep=""), "\n")
}
}
