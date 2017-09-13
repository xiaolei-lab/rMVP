MVP.Version <-
function(){
##############################################################################################
# MVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS
# Designed by Xiaolei Liu
# Writen by Lilin Yin, Zhiwu Zhang, and Xiaolei Liu
# Build date: Aug 30, 2017
# Last update: May 25, 2017
##############################################################################################
version <- "0.1"
authors <- "Lilin Yin, Haohao Zhang, Zhiwu Zhang, Xinyun Li, Shuhong Zhao, and Xiaolei Liu"
cat(paste("#", paste(rep("-", 38), collapse=""), "Welcome to MVP", paste(rep("-", 38), collapse=""), "#",sep=""), "\n")

cat("#   ", "A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS    #","\n")
cat(paste("# ", "Version: ", version, paste(rep(" ", 80-nchar(version)), collapse=""), "#", sep=""), "\n")
cat(paste("#", "Authors:", authors, paste(rep(" ", 78-nchar(authors)), collapse=""),"#", "\n"))
cat(paste("#", "Contact: xiaoleiliu@mail.hzau.edu.cn", paste(rep(" ", 51), collapse=""), "#", "\n"))
cat(paste("#",  paste(rep("-", 90), collapse=""), "#", sep=""), "\n")
}

